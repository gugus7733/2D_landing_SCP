function prob = build_subproblem_2d(X_ref, U_ref, P, t_remaining, retry_level, iter)
%BUILD_SUBPROBLEM_2D Constructs the convex QP for the rocket landing.
%
% This function performs the core work of SCP:
% 1. Linearizes the nonlinear 2D rocket dynamics around the reference
%    trajectory (X_ref, U_ref).
% 2. Formulates the cost function (penalizing fuel, control effort, etc.).
% 3. Sets up equality constraints (dynamics) and inequality constraints
%    (control limits, trust regions).
% 4. Enhanced: Uses time-dependent slack weights and bounds
% 5. Constraint softening: Adds large-penalty slacks for first 2 iterations

% Handle optional parameters for backward compatibility
if nargin < 4
    t_remaining = P.T; % Use full horizon if not specified
end
if nargin < 5
    retry_level = 0; % No retry if not specified
end
if nargin < 6
    iter = 1; % Default to first iteration if not specified
end

N = P.N;
n_x = P.n_states;
n_u = P.n_controls;
dt = P.dt;

% --- Decision Variable Vector 'z' ---
% z = [vec(X); vec(U); vec(s_v); vec(s_w); vec(s_T_upper); vec(s_T_lower); 
%      vec(s_delta_upper); vec(s_delta_lower); vec(s_terminal)]
% where X is states, U is controls, s_v/s_w are dynamics slacks,
% s_T/s_delta are constraint violation slacks, s_terminal is terminal slack
prob.idx.X = 1:(n_x * (N + 1));
prob.idx.U = prob.idx.X(end) + (1:(n_u * N));
prob.idx.s_v = prob.idx.U(end) + (1:(2 * N)); % For vx, vy dynamics
prob.idx.s_w = prob.idx.s_v(end) + (1:(1 * N)); % For omega dynamics

% New constraint violation slacks (for softening during first 2 iterations)
prob.idx.s_T_upper = prob.idx.s_w(end) + (1:N);        % Thrust upper bound violations
prob.idx.s_T_lower = prob.idx.s_T_upper(end) + (1:N);  % Thrust lower bound violations
prob.idx.s_delta_upper = prob.idx.s_T_lower(end) + (1:N); % Gimbal upper bound violations
prob.idx.s_delta_lower = prob.idx.s_delta_upper(end) + (1:N); % Gimbal lower bound violations
prob.idx.s_terminal = prob.idx.s_delta_lower(end) + (1:(n_x-1)); % Terminal constraint violations

n_vars = prob.idx.s_terminal(end);

% --- Cost Function: J = z' H z + f' z ---
H = spalloc(n_vars, n_vars, 4*n_u*N);
f = zeros(n_vars, 1);

% Control magnitude penalties
u_idx = reshape(prob.idx.U, n_u, N);
H(u_idx(1,:), u_idx(1,:)) = H(u_idx(1,:), u_idx(1,:)) + P.w_T * speye(N);
H(u_idx(2,:), u_idx(2,:)) = H(u_idx(2,:), u_idx(2,:)) + P.w_delta * speye(N);

% Control rate penalties (d/dt U)Constraint matrix
D = spdiags([-ones(N-1,1), ones(N-1,1)], [0, 1], N-1, N) / dt;
H_dT = D' * D * P.w_dT;
H_ddelta = D' * D * P.w_ddelta;
H(u_idx(1,:), u_idx(1,:)) = H(u_idx(1,:), u_idx(1,:)) + H_dT;
H(u_idx(2,:), u_idx(2,:)) = H(u_idx(2,:), u_idx(2,:)) + H_ddelta;

% Angular velocity magnitude penalty
if P.w_omega > 0
    x_indices = reshape(prob.idx.X, n_x, N+1);
    omega_idx = x_indices(6, :);  % ω is 6th state at all time points
    H(omega_idx, omega_idx) = H(omega_idx, omega_idx) + P.w_omega * speye(N+1);
end

% Enhanced slack variable penalty (time-dependent L1 norm)
% Use enhanced slack management
w_slack_current = slack_management_utils('compute_slack_weight', t_remaining, P.T, P);
f(prob.idx.s_v) = w_slack_current;
f(prob.idx.s_w) = w_slack_current;

% Large penalty for constraint violation slacks (active only in first 2 iterations)
constraint_slack_penalty = 1e6; % Very large penalty to discourage use
if iter <= 2 || retry_level == 0
    % Apply large penalties to all constraint violation slacks
    f(prob.idx.s_T_upper) = constraint_slack_penalty;
    f(prob.idx.s_T_lower) = constraint_slack_penalty;
    f(prob.idx.s_delta_upper) = constraint_slack_penalty;
    f(prob.idx.s_delta_lower) = constraint_slack_penalty;
    f(prob.idx.s_terminal) = constraint_slack_penalty;
else
    % In later iterations, set infinite penalty (effectively disable slacks)
    f(prob.idx.s_T_upper) = inf;
    f(prob.idx.s_T_lower) = inf;
    f(prob.idx.s_delta_upper) = inf;
    f(prob.idx.s_delta_lower) = inf;
    f(prob.idx.s_terminal) = inf;
end

% --- Dynamics Constraints: Aeq z = beq ---
% Allocate rows for: initial (7) + dynamics (7*N) + terminal (6)
n_eq = n_x * (N + 1) + (n_x - 1);  % 7*(N+1) + 6 = 7*N + 13 total rows
Aeq = spalloc(n_eq, n_vars, n_eq * (n_x + n_u) * 2);
beq = zeros(n_eq, 1);

% fprintf('  Constraint matrix: %d equations, %d variables\n', n_eq, n_vars);
eq_idx = 1:n_x;

% Initial conditions
x_idx_k = 1:n_x;
Aeq(eq_idx, x_idx_k) = speye(n_x);
beq(eq_idx) = P.X0;

% Loop through each time step to build linearized dynamics
for k = 1:N
    % State and control indices for this step (k) and the next (k+1)
    x_idx_k   = (k-1)*n_x + (1:n_x);
    x_idx_kp1 = k*n_x + (1:n_x);
    u_idx_k   = prob.idx.U(1) + (k-1)*n_u + (0:n_u-1);
    sv_idx_k  = prob.idx.s_v(1) + (k-1)*2 + (0:1);
    sw_idx_k  = prob.idx.s_w(1) + (k-1);

    % Get reference state/control and pre-calculate values
    Xk_ref = X_ref(:, k);
    Uk_ref = U_ref(:, k);
    [Ak, Bk, Sk, wk] = get_linearized_dynamics(Xk_ref, Uk_ref, P);

    % Create constraint matrix for this step
    % Form: X_{k+1} = Ak*X_k + Bk*U_k + Sk*S_k + wk
    % --> -Ak*X_k - Bk*U_k + X_{k+1} - Sk*S_k = wk
    eq_idx = n_x*k + (1:n_x);
    Aeq(eq_idx, x_idx_k)   = -Ak;
    Aeq(eq_idx, u_idx_k)   = -Bk;
    Aeq(eq_idx, x_idx_kp1) = speye(n_x);
    Aeq(eq_idx, [sv_idx_k, sw_idx_k]) = -Sk;
    beq(eq_idx) = wk;
end

% Terminal constraints (placed AFTER all dynamics constraints)
target_state = [P.x_target; P.y_target; P.vx_target; P.vy_target; ...
                P.theta_target; P.omega_target; 0]; % mass is free
x_idx_N = (N*n_x) + (1:n_x);  % Final state variables
% Use rows after all dynamics: 7*(N+1) + (1:6)
terminal_eq_rows = (n_x*(N+1)) + (1:n_x-1); % All states except mass
Aeq(terminal_eq_rows, x_idx_N(1:n_x-1)) = speye(n_x-1);

% Add terminal slack variables with -I coefficient (softening during first 2 iterations)
if iter <= 2 || retry_level == 0
    Aeq(terminal_eq_rows, prob.idx.s_terminal) = -speye(n_x-1);
end

beq(terminal_eq_rows) = target_state(1:n_x-1);

% fprintf('  Dynamics rows: 1-%d, Terminal rows: %d-%d\n', ...
%         n_x*(N+1), terminal_eq_rows(1), terminal_eq_rows(end));

% --- Bounds and Inequalities: lb <= z <= ub, Aineq z <= bineq ---
lb = -inf(n_vars, 1);
ub = inf(n_vars, 1);

% State bounds (mass)
x_indices = reshape(prob.idx.X, n_x, N+1);
lb(x_indices(7,:)') = P.m_dry;

% Control bounds (modified for softening in first 2 iterations)
u_indices = reshape(prob.idx.U, n_u, N);

if iter <= 2 || retry_level == 0
    % In first 2 iterations: implement as soft inequality constraints
    % Keep control bounds wide, actual constraints will be in Aineq
    lb(u_indices(1,:)') = -inf;
    ub(u_indices(1,:)') = inf;
    lb(u_indices(2,:)') = -inf;
    ub(u_indices(2,:)') = inf;
else
    % Later iterations: use hard bounds
    lb(u_indices(1,:)') = P.T_min;
    ub(u_indices(1,:)') = P.T_max;
    lb(u_indices(2,:)') = -P.delta_max;
    ub(u_indices(2,:)') = P.delta_max;
end

% Enhanced slack bounds (positive with adaptive upper bounds)
lb(prob.idx.s_v) = 0;
lb(prob.idx.s_w) = 0;

% Positive lower bounds for constraint violation slacks
lb(prob.idx.s_T_upper) = 0;
lb(prob.idx.s_T_lower) = 0;
lb(prob.idx.s_delta_upper) = 0;
lb(prob.idx.s_delta_lower) = 0;
lb(prob.idx.s_terminal) = 0;

% Set adaptive upper bounds for slack variables
if isfield(P, 'slack_max_initial') && isfield(P, 'slack_max_terminal')
    % Use enhanced slack management with retry-aware bounds
    s_max_vector = slack_management_utils('compute_slack_bounds', t_remaining, P.T, P, retry_level);
    
    % Apply bounds to velocity slack variables (s_v: [vx; vy] for each time step)
    s_v_indices = reshape(prob.idx.s_v, 2, N);
    ub(s_v_indices(1,:)') = s_max_vector(1); % vx slack bounds
    ub(s_v_indices(2,:)') = s_max_vector(2); % vy slack bounds
    
    % Apply bounds to angular velocity slack variables (s_w: omega for each time step)
    ub(prob.idx.s_w) = s_max_vector(3); % omega slack bounds
    
    % Log slack bounds if this is a retry
    if retry_level > 0
        slack_management_utils('log_slack_status', t_remaining, P.T, P, 0, retry_level);
    end
end

% Trust region bounds on controls
ub(u_indices(1,:)') = min(ub(u_indices(1,:)'), U_ref(1,:)' + P.trust_T);
lb(u_indices(1,:)') = max(lb(u_indices(1,:)'), U_ref(1,:)' - P.trust_T);
ub(u_indices(2,:)') = min(ub(u_indices(2,:)'), U_ref(2,:)' + P.trust_delta);
lb(u_indices(2,:)') = max(lb(u_indices(2,:)'), U_ref(2,:)' - P.trust_delta);

% Trust region bounds on states (matching 1D model approach)
if isfield(P, 'trust_vx') && isfield(P, 'trust_vy') && isfield(P, 'trust_omega')
    % Velocity trust regions
    ub(x_indices(3,:)') = min(ub(x_indices(3,:)'), X_ref(3,:)' + P.trust_vx);
    lb(x_indices(3,:)') = max(lb(x_indices(3,:)'), X_ref(3,:)' - P.trust_vx);
    ub(x_indices(4,:)') = min(ub(x_indices(4,:)'), X_ref(4,:)' + P.trust_vy);
    lb(x_indices(4,:)') = max(lb(x_indices(4,:)'), X_ref(4,:)' - P.trust_vy);
    % Angular velocity trust regions
    ub(x_indices(6,:)') = min(ub(x_indices(6,:)'), X_ref(6,:)' + P.trust_omega);
    lb(x_indices(6,:)') = max(lb(x_indices(6,:)'), X_ref(6,:)' - P.trust_omega);
end

% --- Inequality Constraints for Soft Control Bounds ---
if iter <= 2 || retry_level == 0
    % Create inequality constraints: A_ineq * z <= b_ineq
    % For thrust: T_k <= T_max + s_T_upper_k  =>  T_k - s_T_upper_k <= T_max
    %             T_k >= T_min - s_T_lower_k  =>  -T_k - s_T_lower_k <= -T_min
    % For gimbal: delta_k <= delta_max + s_delta_upper_k  =>  delta_k - s_delta_upper_k <= delta_max
    %             delta_k >= -delta_max - s_delta_lower_k  =>  -delta_k - s_delta_lower_k <= delta_max
    
    n_ineq = 4 * N; % 4 inequalities per time step (T_upper, T_lower, delta_upper, delta_lower)
    Aineq = spalloc(n_ineq, n_vars, 2 * n_ineq);
    bineq = zeros(n_ineq, 1);
    
    % Build inequality constraints for each time step
    for k = 1:N
        u_idx_k = prob.idx.U(1) + (k-1)*n_u + (0:n_u-1);
        
        % Row indices for this time step's constraints
        row_base = (k-1)*4;
        
        % Thrust upper bound: T_k - s_T_upper_k <= T_max
        Aineq(row_base + 1, u_idx_k(1)) = 1;
        Aineq(row_base + 1, prob.idx.s_T_upper(k)) = -1;
        bineq(row_base + 1) = P.T_max;
        
        % Thrust lower bound: -T_k - s_T_lower_k <= -T_min
        Aineq(row_base + 2, u_idx_k(1)) = -1;
        Aineq(row_base + 2, prob.idx.s_T_lower(k)) = -1;
        bineq(row_base + 2) = -P.T_min;
        
        % Gimbal upper bound: delta_k - s_delta_upper_k <= delta_max
        Aineq(row_base + 3, u_idx_k(2)) = 1;
        Aineq(row_base + 3, prob.idx.s_delta_upper(k)) = -1;
        bineq(row_base + 3) = P.delta_max;
        
        % Gimbal lower bound: -delta_k - s_delta_lower_k <= delta_max
        Aineq(row_base + 4, u_idx_k(2)) = -1;
        Aineq(row_base + 4, prob.idx.s_delta_lower(k)) = -1;
        bineq(row_base + 4) = P.delta_max;
    end
else
    % Later iterations: no inequality constraints (use hard bounds)
    Aineq = [];
    bineq = [];
end

prob.H = H; prob.f = f;
prob.Aeq = Aeq; prob.beq = beq;
prob.Aineq = Aineq; prob.bineq = bineq;
prob.lb = lb; prob.ub = ub;
prob.P = P;
end


function [Ad,Bd,Sd,wd] = get_linearized_dynamics(Xk,Uk,P)
%GET_LINEARIZED_DYNAMICS_CASADI  Linear‐affine discretised model via CasADi AD
%
%   Continuous model (around reference (Xk,Uk)):
%       ẋ = f(x,u)  ≈  Ac·(x-x̂) + Bc·(u-û) + wc
%
%   Discretised with forward Euler:
%       x⁺ = Ad x + Bd u + Sd s + wd
%
%   Returns:
%       Ad, Bd : discrete Jacobians
%       Sd     : slack gain matrix (unchanged from hand version)
%       wd     : discrete affine term

% ---------- 1.  Initialise CasADi symbols once (persistent) ----------
import casadi.*

persistent f_sym A_sym B_sym W_sym stateDim ctrlDim
if isempty(f_sym)
    % Dimensions -------------------------------------------------------
    stateDim = P.n_states;   % 7 for [x y vx vy th om m]
    ctrlDim  = P.n_controls; % 2 for [T delta]

    % Symbolic variables ----------------------------------------------
    X  = SX.sym('X', stateDim);      % [x y vx vy th om m]
    U  = SX.sym('U', ctrlDim);       % [T delta]

    % Unpack for readability ------------------------------------------
    x  = X(1);  y  = X(2);
    vx = X(3);  vy = X(4);
    th = X(5);  om = X(6);
    m  = X(7);

    T     = U(1);
    delta = U(2);

    % ---------- 2.  Dynamics definition ------------------------------
    cth = cos(th); sth = sin(th);
    cdel = cos(delta); sdel = sin(delta);

    % Kinematics
    x_dot  = vx;
    y_dot  = vy;
    th_dot = om;

    % Atmosphere (standard ISA)
    [~,~,~,rho] = atmosisa(Xk(2));
    V_sq  = vx^2 + vy^2;
    q     = 0.5*rho*P.A_ref*V_sq;

    % Angle of attack (small-angle approx; CasADi differentiable)
    v_I      = [-vx; -vy];                  % into the flow
    T_B_I    = [ cth, sth; -sth, cth];
    v_aero_B = T_B_I * v_I;
    alpha    = atan2(v_aero_B(1), v_aero_B(2));

    % Aero forces in body & inertial frames
    F_ax_B =  q * P.Cd;                  % along body y-axis
    F_n_B  = -q * P.C_N_alpha * alpha;   % along body x-axis
    T_I_B  = [ cth, -sth; sth, cth];
    F_aero_I = T_I_B * [F_n_B; F_ax_B];

    % Thrust (body → gimballed → inertial)
    TBthrust  = [cos(delta), -sin(delta);
                 sin(delta),  cos(delta)];
    F_thrust_I = T_I_B * (TBthrust * [0; T]);

    % Translational dynamics
    F_grav  = [0; -P.g0*m];
    v_dot   = (F_thrust_I + F_aero_I + F_grav) / m;

    % Rotational dynamics
    Iyy      = P.Iyy_func(m);
    M_thrust = -T * P.L_com_from_base * sin(delta);
    M_aero   = F_n_B * (P.L_cop_from_base - P.L_com_from_base);
    M_damp   = -P.C_damp * om;
    om_dot   = (M_thrust + M_aero + M_damp) / Iyy;

    % Mass flow
    m_dot = -T / (P.Isp * P.g0);

    % Assemble f(x,u)
    f_cont = [x_dot;
              y_dot;
              v_dot(1);
              v_dot(2);
              th_dot;
              om_dot;
              m_dot];

    % ---------- 3.  Exact Jacobians via AD ----------------------------
    A_sym = Function('A_sym',{X,U},{jacobian(f_cont,X)});
    B_sym = Function('B_sym',{X,U},{jacobian(f_cont,U)});
    W_sym = Function('W_sym',{X,U},{f_cont});   % value of f itself
end

% ---------- 4.  Evaluate at reference point --------------------------
Ac = full(A_sym(Xk,Uk));
Bc = full(B_sym(Xk,Uk));
fc = full(W_sym(Xk,Uk));

% Affine drift term  wc  so that  f(x̂,û) = Ac*x̂ + Bc*û + wc
wc = fc - Ac*Xk - Bc*Uk;

% ---------- 5.  Slack matrix (unchanged) -----------------------------
Sc      = zeros(stateDim,3);
m       = Xk(7);
Iyy     = P.Iyy_func(m);
Sc(3,1) = 1/m;
Sc(4,2) = 1/m;
Sc(6,3) = 1/Iyy;

% ---------- 6.  Discretisation (forward Euler) -----------------------
Ad = eye(stateDim) + Ac * P.dt;
Bd =            Bc * P.dt;
Sd =            Sc * P.dt;
wd =            wc * P.dt;
end
