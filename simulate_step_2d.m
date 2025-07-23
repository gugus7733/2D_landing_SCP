function new_state = simulate_step_2d(state, control, dt, P)
%SIMULATE_STEP_2D Propagates the 2D rocket dynamics for one time step.
%
% This function implements the true nonlinear equations of motion.

% Unpack state and control
% State unpacking: [x; y; vx; vy; theta; omega; m]
x=state(1); y=state(2); vx=state(3); vy=state(4); th=state(5); om=state(6); m=state(7);
T_cmd=control(1); delta_cmd=control(2);

% --- Enforce physical limits ---
T = min(max(T_cmd, P.T_min), P.T_max);
delta = min(max(delta_cmd, -P.delta_max), P.delta_max);
if m <= P.m_dry, T = 0; end % Flameout

% --- Mass flow rate ---
mdot = T / (P.Isp * P.g0);

% --- Forces and Moments (in Inertial Frame) ---
% Gravity (always in inertial frame, pointing down)
F_grav = [0; -m * P.g0];

% Thrust (computed via proper frame transformations)
% 1. Define thrust vector in body frame along yb-axis
thrust_vector_B = [0; T];
% 2. Apply gimbal rotation in body frame
T_B_thrust = [cos(delta), -sin(delta); sin(delta), cos(delta)];
thrust_gimbaled_B = T_B_thrust * thrust_vector_B;
% 3. Transform to inertial frame
T_I_B = [cos(th), -sin(th); sin(th), cos(th)];
F_thrust = T_I_B * thrust_gimbaled_B;


% Aerodynamics
[~, ~, ~, rho] = atmosisa(y); % Use altitude y for standard atmosphere
V_sq = vx^2 + vy^2;
q = 0.5 * rho * P.A_ref * V_sq;

if V_sq > 1e-3
    % Velocity vector in inertial frame
    v_I = [vx; vy];
    % Aerodynamic forces are opposed to velocity, normalized
    v_aero_I = -v_I / sqrt(V_sq);
    % Transform to body frame to compute angle of attack
    T_B_I = [cos(th), sin(th); -sin(th), cos(th)];
    v_aero_B = T_B_I * v_aero_I;
    alpha = atan2(v_aero_B(1), v_aero_B(2));
    
    % Compute aerodynamic forces in body frame
    F_axial_B = q * P.Cd;      % Along yb-axis (drag)
    F_normal_B = -q * P.C_N_alpha * alpha;  % Along xb-axis (side force)
    F_aero_B = [F_normal_B; F_axial_B];
    
    % Transform to inertial frame
    F_aero = T_I_B * F_aero_B;
else
    % Velocity vector in inertial frame
    v_I = [0; 0];
    % Aerodynamic forces are opposed to velocity, normalized
    v_aero_I = v_I;
    v_aero_B = v_aero_I;
    alpha = 0;
    
    % Compute aerodynamic forces in body frame
    F_axial_B = 0;      % Along yb-axis (drag)
    F_normal_B = 0;  % Along xb-axis (side force)
    F_aero_B = [F_normal_B; F_axial_B];
    F_aero = 0;
end

% --- Dynamics ---
% Translational
F_total = F_thrust + F_aero + F_grav;
accel = F_total / m;

% Rotational
Iyy = P.Iyy_func(m);
% Thrust moment: T * L_com * sin(delta) but with correct sign convention
% Positive delta creates negative moment (nose down tendency)
M_thrust = -T * P.L_com_from_base * sin(delta);
% Aerodynamic moment: normal force * moment arm
% Positive alpha creates stabilizing moment (negative)
M_aero = F_aero_B(1) * (P.L_cop_from_base - P.L_com_from_base);

% Aerodynamic damping moment (proportional to omega, applied at CoP)
if ~isfield(P, 'C_damp')
    P.C_damp = 0.0; % Default: no damping if not set
end
M_damp = -P.C_damp * om; % Damping opposes rotation

omega_dot = (M_thrust + M_aero + M_damp) / Iyy;

% --- Integration (Simple Euler method) ---
% Integration (Simple Euler method)
x_new     = x + vx * dt;
y_new     = y + vy * dt;
vx_new    = vx + accel(1) * dt;
vy_new    = vy + accel(2) * dt;
th_new    = th + om * dt;
om_new    = om + omega_dot * dt;
m_new     = m - mdot * dt;

new_state = [x_new; y_new; vx_new; vy_new; th_new; om_new; max(m_new, P.m_dry)];
end