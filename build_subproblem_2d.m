function prob = build_subproblem_2d(X_ref, U_ref, P)
%BUILD_SUBPROBLEM_2D Constructs the convex QP for the rocket landing.
%
% This function performs the core work of SCP:
% 1. Linearizes the nonlinear 2D rocket dynamics around the reference
%    trajectory (X_ref, U_ref).
% 2. Formulates the cost function (penalizing fuel, control effort, etc.).
% 3. Sets up equality constraints (dynamics) and inequality constraints
%    (control limits, trust regions).

N = P.N;
n_x = P.n_states;
n_u = P.n_controls;
dt = P.dt;

% --- Decision Variable Vector 'z' ---
% z = [vec(X); vec(U); vec(s_v); vec(s_w)]
% where X is states, U is controls, s_v/s_w are slacks
prob.idx.X = 1:(n_x * (N + 1));
prob.idx.U = prob.idx.X(end) + (1:(n_u * N));
prob.idx.s_v = prob.idx.U(end) + (1:(2 * N)); % For vx, vy dynamics
prob.idx.s_w = prob.idx.s_v(end) + (1:(1 * N)); % For omega dynamics
n_vars = prob.idx.s_w(end);

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

% Slack variable penalty (L1 norm)
f(prob.idx.s_v) = P.w_slack;
f(prob.idx.s_w) = P.w_slack;

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
beq(terminal_eq_rows) = target_state(1:n_x-1);

% fprintf('  Dynamics rows: 1-%d, Terminal rows: %d-%d\n', ...
%         n_x*(N+1), terminal_eq_rows(1), terminal_eq_rows(end));

% --- Bounds and Inequalities: lb <= z <= ub, Aineq z <= bineq ---
lb = -inf(n_vars, 1);
ub = inf(n_vars, 1);

% State bounds (mass)
x_indices = reshape(prob.idx.X, n_x, N+1);
lb(x_indices(7,:)') = P.m_dry;

% Control bounds
u_indices = reshape(prob.idx.U, n_u, N);
lb(u_indices(1,:)') = P.T_min;
ub(u_indices(1,:)') = P.T_max;
lb(u_indices(2,:)') = -P.delta_max;
ub(u_indices(2,:)') = P.delta_max;

% Slack bounds (must be positive)
lb(prob.idx.s_v) = 0;
lb(prob.idx.s_w) = 0;

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

prob.H = H; prob.f = f;
prob.Aeq = Aeq; prob.beq = beq;
prob.Aineq = []; prob.bineq = [];
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
