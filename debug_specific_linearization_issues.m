function debug_specific_linearization_issues(X_ref, U_ref, initial_state, T_horizon, N, P)
%DEBUG_SPECIFIC_LINEARIZATION_ISSUES Tests specific suspected issues in linearization
%
% This function focuses on specific potential problems identified in the 
% get_linearized_dynamics function

fprintf('\n=== SPECIFIC LINEARIZATION ISSUE DEBUGGING ===\n');

% Test point
X_test = X_ref(:, 1);
U_test = U_ref(:, 1);

P_scp = P;
P_scp.X0 = initial_state;
P_scp.T = T_horizon;
P_scp.N = N;
P_scp.dt = T_horizon / N;
P_scp.max_iters = P.max_iters_scp;

% State and Control dimensions
P_scp.n_states = 7;
P_scp.n_controls = 2;
P = P_scp;

% Extract values
x=X_test(1); y=X_test(2); vx=X_test(3); vy=X_test(4); 
th=X_test(5); om=X_test(6); m=X_test(7);
T=U_test(1); delta=U_test(2);

fprintf('Test point: x=%.1f, y=%.1f, vx=%.1f, vy=%.1f, th=%.2f°, m=%.0f kg\n', ...
    x, y, vx, vy, rad2deg(th), m);
fprintf('Controls: T=%.1f kN, delta=%.2f°\n', T/1e3, rad2deg(delta));

%% Issue 1: Check affine term calculation
fprintf('\n1. AFFINE TERM CALCULATION:\n');

% Get linearized dynamics
[Ad, Bd, Sd, wd] = get_linearized_dynamics(X_test, U_test, P);

% The continuous-time affine term should be: f(x_ref, u_ref) - A*x_ref - B*u_ref
% where f is the nonlinear dynamics RHS

% Compute nonlinear dynamics at reference point
f_nonlinear = compute_dynamics_rhs(X_test, U_test, P);

% Extract continuous-time matrices (undo discretization)
Ac = (Ad - eye(7)) / P.dt;
Bc = Bd / P.dt;
wc = wd / P.dt;

% What the affine term should be
wc_correct = f_nonlinear - Ac * X_test - Bc * U_test;

% Compare
affine_error = norm(wc - wc_correct);
fprintf('  Computed affine term norm: %.6f\n', norm(wc));
fprintf('  Correct affine term norm:  %.6f\n', norm(wc_correct));
fprintf('  Affine term error:         %.2e\n', affine_error);

if affine_error > 1e-10
    fprintf('  *** AFFINE TERM ERROR DETECTED ***\n');
    fprintf('  Difference breakdown:\n');
    state_names = {'x', 'y', 'vx', 'vy', 'theta', 'omega', 'm'};
    for i = 1:7
        fprintf('    %s: computed=%.6f, correct=%.6f, error=%.2e\n', ...
            state_names{i}, wc(i), wc_correct(i), abs(wc(i) - wc_correct(i)));
    end
end

%% Issue 2: Check thrust force derivatives
fprintf('\n2. THRUST FORCE DERIVATIVES:\n');

% Manual calculation of thrust force in inertial frame
cth = cos(th); sth = sin(th);
cdelta = cos(delta); sdelta = sin(delta);

% Thrust vector transformation: 
% 1. Start with [0; T] in body frame
% 2. Apply gimbal: [T*sin(delta); T*cos(delta)]
% 3. Transform to inertial: T_I_B * [T*sin(delta); T*cos(delta)]
thrust_I_manual = [T*(cth*sdelta - sth*cdelta); T*(sth*sdelta + cth*cdelta)];

% Compare with simulate_step_2d calculation
thrust_vector_B = [0; T];
T_B_thrust = [cos(delta), -sin(delta); sin(delta), cos(delta)];
thrust_gimbaled_B = T_B_thrust * thrust_vector_B;
T_I_B = [cos(th), -sin(th); sin(th), cos(th)];
thrust_I_simulate = T_I_B * thrust_gimbaled_B;

thrust_calc_error = norm(thrust_I_manual - thrust_I_simulate);
fprintf('  Manual calculation:    [%.1f, %.1f] N\n', thrust_I_manual);
fprintf('  Simulate calculation:  [%.1f, %.1f] N\n', thrust_I_simulate);
fprintf('  Calculation error:     %.2e\n', thrust_calc_error);

% Check partial derivatives
eps = 1e-8;

% d(F_thrust)/dT
T_plus = T + eps;
thrust_I_plus = [T_plus*(cth*sdelta - sth*cdelta); T_plus*(sth*sdelta + cth*cdelta)];
dF_dT_numerical = (thrust_I_plus - thrust_I_manual) / eps;
dF_dT_analytical = [(cth*sdelta - sth*cdelta); (sth*sdelta + cth*cdelta)];

dF_dT_error = norm(dF_dT_numerical - dF_dT_analytical);
fprintf('  dF/dT numerical:   [%.6f, %.6f]\n', dF_dT_numerical);
fprintf('  dF/dT analytical:  [%.6f, %.6f]\n', dF_dT_analytical);
fprintf('  dF/dT error:       %.2e\n', dF_dT_error);

% d(F_thrust)/ddelta  
delta_plus = delta + eps;
cplus = cos(delta_plus); splus = sin(delta_plus);
thrust_I_delta_plus = [T*(cth*splus - sth*cplus); T*(sth*splus + cth*cplus)];
dF_ddelta_numerical = (thrust_I_delta_plus - thrust_I_manual) / eps;
dF_ddelta_analytical = T * [(cth*cdelta + sth*sdelta); (sth*cdelta - cth*sdelta)];

dF_ddelta_error = norm(dF_ddelta_numerical - dF_ddelta_analytical);
fprintf('  dF/ddelta numerical:   [%.6f, %.6f]\n', dF_ddelta_numerical);
fprintf('  dF/ddelta analytical:  [%.6f, %.6f]\n', dF_ddelta_analytical);
fprintf('  dF/ddelta error:       %.2e\n', dF_ddelta_error);

%% Issue 3: Check moment derivatives
fprintf('\n3. MOMENT DERIVATIVES:\n');

Iyy = P.Iyy_func(m);

% Moment from thrust: M = -T * L_com * sin(delta)
M_thrust = -T * P.L_com_from_base * sin(delta);

% Check derivatives
dM_dT_analytical = -P.L_com_from_base * sin(delta);
dM_ddelta_analytical = -T * P.L_com_from_base * cos(delta);

% Numerical check
T_plus = T + eps;
M_thrust_T_plus = -T_plus * P.L_com_from_base * sin(delta);
dM_dT_numerical = (M_thrust_T_plus - M_thrust) / eps;

delta_plus = delta + eps;
M_thrust_delta_plus = -T * P.L_com_from_base * sin(delta_plus);
dM_ddelta_numerical = (M_thrust_delta_plus - M_thrust) / eps;

fprintf('  dM/dT numerical:   %.6f\n', dM_dT_numerical);
fprintf('  dM/dT analytical:  %.6f\n', dM_dT_analytical);
fprintf('  dM/dT error:       %.2e\n', abs(dM_dT_numerical - dM_dT_analytical));

fprintf('  dM/ddelta numerical:   %.6f\n', dM_ddelta_numerical);
fprintf('  dM/ddelta analytical:  %.6f\n', dM_ddelta_analytical);
fprintf('  dM/ddelta error:       %.2e\n', abs(dM_ddelta_numerical - dM_ddelta_analytical));

% Check what's in the B matrix
fprintf('  B matrix moment entries: Bc(6,1)=%.6f, Bc(6,2)=%.6f\n', ...
    Bc(6,1)*P.dt, Bc(6,2)*P.dt);  % Convert back to continuous time
fprintf('  Expected: %.6f, %.6f\n', dM_dT_analytical/Iyy, dM_ddelta_analytical/Iyy);

%% Issue 4: Discretization accuracy
fprintf('\n4. DISCRETIZATION ACCURACY:\n');

dt_values = [0.01, 0.05, 0.1, 0.2];
for dt_test = dt_values
    % Forward Euler: X_{k+1} = X_k + dt * f(X_k, U_k)
    X_euler = X_test + dt_test * f_nonlinear;
    
    % Higher order (RK4) for comparison
    k1 = f_nonlinear;
    k2 = compute_dynamics_rhs(X_test + dt_test*k1/2, U_test, P);
    k3 = compute_dynamics_rhs(X_test + dt_test*k2/2, U_test, P);
    k4 = compute_dynamics_rhs(X_test + dt_test*k3, U_test, P);
    X_rk4 = X_test + dt_test * (k1 + 2*k2 + 2*k3 + k4) / 6;
    
    % True nonlinear simulation
    X_true = simulate_step_2d(X_test, U_test, dt_test, P);
    
    euler_error = norm(X_euler - X_true);
    rk4_error = norm(X_rk4 - X_true);
    
    fprintf('  dt=%.3f: Euler error=%.2e, RK4 error=%.2e\n', dt_test, euler_error, rk4_error);
end

fprintf('\n=== SUMMARY OF SPECIFIC ISSUES ===\n');
if affine_error > 1e-10
    fprintf('*** CRITICAL: Affine term calculation error ***\n');
end
if dF_dT_error > 1e-10 || dF_ddelta_error > 1e-10
    fprintf('*** CRITICAL: Thrust force derivative errors ***\n');
end
if abs(dM_dT_numerical - dM_dT_analytical) > 1e-10 || abs(dM_ddelta_numerical - dM_ddelta_analytical) > 1e-10
    fprintf('*** CRITICAL: Moment derivative errors ***\n');
end

end

function f = compute_dynamics_rhs(X, U, P)
%COMPUTE_DYNAMICS_RHS Same as in validate_linearization.m

% Extract state
x=X(1); y=X(2); vx=X(3); vy=X(4); th=X(5); om=X(6); m=X(7);
T=U(1); delta=U(2);

% Apply limits
T = min(max(T, P.T_min), P.T_max);
delta = min(max(delta, -P.delta_max), P.delta_max);
if m <= P.m_dry, T = 0; end

% Mass flow rate
mdot = T / (P.Isp * P.g0);

% Forces computation (same as simulate_step_2d but continuous-time)
F_grav = [0; -m * P.g0];

% Thrust forces
thrust_vector_B = [0; T];
T_B_thrust = [cos(delta), -sin(delta); sin(delta), cos(delta)];
thrust_gimbaled_B = T_B_thrust * thrust_vector_B;
T_I_B = [cos(th), -sin(th); sin(th), cos(th)];
F_thrust = T_I_B * thrust_gimbaled_B;

% Aerodynamics
[~, ~, ~, rho] = atmosisa(y);
V_sq = vx^2 + vy^2;
q = 0.5 * rho * P.A_ref * V_sq;

if V_sq > 1e-3
    v_I = [vx; vy];
    v_aero_I = -v_I / sqrt(V_sq);
    T_B_I = [cos(th), sin(th); -sin(th), cos(th)];
    v_aero_B = T_B_I * v_aero_I;
    alpha = atan2(v_aero_B(1), v_aero_B(2));
    
    F_axial_B = q * P.Cd;
    F_normal_B = -q * P.C_N_alpha * alpha;
    F_aero_B = [F_normal_B; F_axial_B];
    F_aero = T_I_B * F_aero_B;
else
    alpha = 0;
    F_aero = [0; 0];
end

% Translational dynamics
F_total = F_thrust + F_aero + F_grav;
accel = F_total / m;

% Rotational dynamics
Iyy = P.Iyy_func(m);
M_thrust = -T * P.L_com_from_base * sin(delta);
M_aero = F_aero_B(1) * (P.L_cop_from_base - P.L_com_from_base);
if ~isfield(P, 'C_damp'), P.C_damp = 0.0; end
M_damp = -P.C_damp * om;
omega_dot = (M_thrust + M_aero + M_damp) / Iyy;

% Return state derivatives
f = [vx; vy; accel(1); accel(2); om; omega_dot; -mdot];

end

function [Ad, Bd, Sd, wd] = get_linearized_dynamics(Xk, Uk, P)
% GET_LINEARIZED_DYNAMICS
% Computes the continuous-time linearized dynamics matrices (A, B) and
% the affine term (w), then discretizes them.
%
% X_dot = A*X + B*U + S*s + w

n_x = P.n_states;
n_u = P.n_controls;

% Extract reference state and control
x=Xk(1); y=Xk(2); vx=Xk(3); vy=Xk(4); th=Xk(5); om=Xk(6); m=Xk(7);
T=Uk(1); delta=Uk(2);

% --- Continuous-time dynamics linearization: x_dot = f(x,u) ---
% A = df/dx, B = df/du, w = f(x_ref, u_ref) - A*x_ref - B*u_ref
Ac = zeros(n_x, n_x);
Bc = zeros(n_x, n_u);

% Kinematic equations (already linear)
Ac(1,3) = 1; % x_dot = vx
Ac(2,4) = 1; % y_dot = vy
Ac(5,6) = 1; % th_dot = omega

% Mass dynamics
Bc(7,1) = -1 / (P.Isp * P.g0); % dm/dT

    % Force and Moment calculations (at reference point)
    [~, ~, ~, rho] = atmosisa(y); % Use altitude y for standard atmosphere
    V_sq = vx^2 + vy^2;
    q = 0.5 * rho * P.A_ref * V_sq;

% Compute angle of attack using proper frame transformations
if V_sq > 1e-3
    % Velocity vector in inertial frame
    v_I = [vx; vy];
    % Aerodynamic forces are opposed to velocity, normalized
    v_aero_I = -v_I / sqrt(V_sq);
    % Transform to body frame to compute angle of attack
    T_B_I = [cos(th), sin(th); -sin(th), cos(th)];
    v_aero_B = T_B_I * v_aero_I;
    alpha = atan2(v_aero_B(1), v_aero_B(2));
else
    alpha = 0;
end

% Aerodynamic forces in body frame
F_ax_B = q * P.Cd;      % Axial drag (along yb-axis)
F_n_B  = -q * P.C_N_alpha * alpha;  % Normal force (along xb-axis)

% Transform aerodynamic forces to inertial frame
T_I_B = [cos(th), -sin(th); sin(th), cos(th)];
F_aero_B = [F_n_B; F_ax_B];
F_aero_I = T_I_B * F_aero_B;

% Thrust forces using proper frame transformations
% Thrust vector in body frame: [0; T]
% Apply gimbal rotation: T_B_thrust = [cos(delta), -sin(delta); sin(delta), cos(delta)]
% Transform to inertial frame: T_I_B
thrust_vector_B = [0; T];
T_B_thrust = [cos(delta), -sin(delta); sin(delta), cos(delta)];
thrust_gimbaled_B = T_B_thrust * thrust_vector_B;
F_thrust_I = T_I_B * thrust_gimbaled_B;

% Rotational dynamics: Iyy * omega_dot = M_thrust + M_aero + M_damp
Iyy = P.Iyy_func(m);
M_thrust = -T * P.L_com_from_base * sin(delta);
M_aero = -q * P.C_N_alpha * alpha * (P.L_cop_from_base - P.L_com_from_base);
if ~isfield(P, 'C_damp')
    P.C_damp = 0.0;
end
M_damp = -P.C_damp * om;
omega_dot = (M_thrust + M_aero + M_damp) / Iyy;

% Translational dynamics: m * v_dot = F_thrust + F_aero + F_gravity
F_thrust_x = F_thrust_I(1);
F_thrust_y = F_thrust_I(2);
F_aero_x = F_aero_I(1);
F_aero_y = F_aero_I(2);
vx_dot = (F_thrust_x + F_aero_x) / m;
vy_dot = (F_thrust_y + F_aero_y) / m - P.g0;

% --- Partial Derivatives (building Ac and Bc) ---
% This is the core of the linearization.
% Partial derivatives of accelerations w.r.t. state and control.
% Note: many derivatives are complex; we use simplified but effective ones.

% d(accel)/d(state) -> Ac
Ac(3,7) = -vx_dot / m; % d(vx_dot)/dm
Ac(4,7) = -(vy_dot + P.g0) / m; % d(vy_dot)/dm
Ac(6,7) = -omega_dot / Iyy * (P.Iyy_func(m+1)-P.Iyy_func(m-1))/2; % d(om_dot)/dm

% Simplification: Aerodynamic cross-coupling terms are often small in this
% regime and are ignored here to maintain robustness of the QP.
% A full analytical derivative would include d(F_aero)/d(vx,vz,th).

% d(accel)/d(control) -> Bc
% For thrust forces: F_thrust = T_I_B * T_B_thrust * [0; T]
cth = cos(th); sth = sin(th);
cdelta = cos(delta); sdelta = sin(delta);

% Thrust force derivatives
% F_thrust_I = T_I_B * [T*sin(delta); T*cos(delta)]
% F_thrust_I = [T*(cth*sdelta - sth*cdelta); T*(sth*sdelta + cth*cdelta)]
Bc(3,1) = (cth*sdelta - sth*cdelta)/m; % d(vx_dot)/dT
Bc(4,1) = (sth*sdelta + cth*cdelta)/m; % d(vy_dot)/dT
Bc(3,2) = T*(cth*cdelta + sth*sdelta)/m; % d(vx_dot)/ddelta
Bc(4,2) = T*(sth*cdelta - cth*sdelta)/m; % d(vy_dot)/ddelta

% Moment derivatives
Bc(6,1) = -P.L_com_from_base*sdelta/Iyy; % d(om_dot)/dT
Bc(6,2) = -T*P.L_com_from_base*cdelta/Iyy; % d(om_dot)/ddelta

% Affine term (or 'offset')
wc = [vx; vy; vx_dot; vy_dot; om; omega_dot; Bc(7,1)*T] - Ac*Xk - Bc*Uk;

% Slack matrix (for virtual controls on accelerations)
Sc = zeros(n_x, 3);
Sc(3,1) = 1/m;   % Slack on vx_dot
Sc(4,2) = 1/m;   % Slack on vy_dot
Sc(6,3) = 1/Iyy; % Slack on omega_dot

% --- Discretization (Forward Euler) ---
Ad = eye(n_x) + Ac * P.dt;
Bd = Bc * P.dt;
Sd = Sc * P.dt;
wd = wc * P.dt;
end