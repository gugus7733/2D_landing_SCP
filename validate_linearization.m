function validate_linearization(X_ref, U_ref, initial_state, T_horizon, N, P)
%VALIDATE_LINEARIZATION Comprehensive validation of SCP linearization accuracy
%
% This function performs detailed validation of the linearized dynamics by:
% 1. Comparing linearized prediction vs true nonlinear dynamics
% 2. Numerical derivative validation of partial derivatives  
% 3. Coordinate transformation verification
% 4. Step size sensitivity analysis
%
% Usage: validate_linearization(X_ref, U_ref, P)
% where X_ref, U_ref are reference trajectories from SCP

fprintf('\n=== LINEARIZATION VALIDATION ===\n');

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

% Test parameters
dt_test = P.dt_scp;
eps = 1e-6;  % For numerical derivatives
% n_test_points = min(5, size(X_ref, 2));  % Test first few points
n_test_points = size(X_ref, 2);  % Test first few points

% Storage for results
validation_results = struct();
validation_results.max_error = 0;
validation_results.avg_error = 0;
validation_results.problematic_states = [];

%% Test 1: Linearized vs Nonlinear Dynamics Comparison
fprintf('\n1. LINEARIZED vs NONLINEAR DYNAMICS:\n');

for k = 1:n_test_points
    if k > size(U_ref, 2), break; end
    
    % Reference point
    X_k = X_ref(:, k);
    U_k = U_ref(:, k);
    
    % Get linearized dynamics matrices
    [A_k, B_k, S_k, w_k] = get_linearized_dynamics(X_k, U_k, P);
    
    % Predict next state using linearized model
    X_pred_linear = A_k * X_k + B_k * U_k + w_k;
    
    % Compute next state using true nonlinear dynamics
    X_pred_nonlinear = simulate_step_2d(X_k, U_k, dt_test, P);
    
    % Compute error
    error = norm(X_pred_linear - X_pred_nonlinear);
    rel_error = error / max(norm(X_pred_nonlinear), 1e-12);
    
    fprintf('  Point %d: ||error|| = %.2e, rel_error = %.2e\n', k, error, rel_error);
    
    % Detailed state-by-state comparison
    state_names = {'x', 'y', 'vx', 'vy', 'theta', 'omega', 'm'};
    for i = 1:length(state_names)
        state_error = abs(X_pred_linear(i) - X_pred_nonlinear(i));
        if state_error > 1e-3  % Significant error threshold
            fprintf('    WARNING: %s error = %.3e (linear=%.3f, nonlinear=%.3f)\n', ...
                state_names{i}, state_error, X_pred_linear(i), X_pred_nonlinear(i));
        end
    end
    
    validation_results.max_error = max(validation_results.max_error, error);
    validation_results.avg_error = validation_results.avg_error + error / n_test_points;
    
    if rel_error > 0.1  % 10% relative error is problematic
        validation_results.problematic_states(end+1) = k;
    end
end

%% Test 1.5: Comprehensive Linear vs Nonlinear Dynamics Visualization
fprintf('\n1.5. COMPREHENSIVE DYNAMICS COMPARISON VISUALIZATION:\n');

% Extended simulation to compare linear vs nonlinear over multiple steps
% n_sim_steps = min(20, size(X_ref, 2)-1);
n_sim_steps = size(X_ref, 2)-1;
if n_sim_steps < 5
    fprintf('   Insufficient reference trajectory for comprehensive comparison.\n');
else
    fprintf('   Simulating %d steps for comprehensive comparison...\n', n_sim_steps);
    
    % Storage for comparison results
    t_sim = (0:n_sim_steps) * dt_test;
    X_linear_sim = zeros(7, n_sim_steps+1);
    X_nonlinear_sim = zeros(7, n_sim_steps+1);
    errors_sim = zeros(1, n_sim_steps);
    
    % Initial condition
    X_linear_sim(:,1) = X_ref(:,1);
    X_nonlinear_sim(:,1) = X_ref(:,1);
    
    % Simulate forward using PROPER validation methodology
    for k = 1:n_sim_steps
        if k > size(U_ref, 2), break; end
        
        % Current states and control
        X_linear_k = X_linear_sim(:, k);
        X_nonlinear_k = X_nonlinear_sim(:, k);
        U_ref_k = U_ref(:, k);
        
        % METHOD 1: Linearize around reference (what SCP actually does)
        X_ref_k = X_ref(:, k);
        [A_ref, B_ref, ~, w_ref] = get_linearized_dynamics(X_ref_k, U_ref_k, P);
        X_linear_sim(:, k+1) = A_ref * X_linear_k + B_ref * U_ref_k + w_ref;
        
        % METHOD 2: True nonlinear prediction
        X_nonlinear_sim(:, k+1) = simulate_step_2d(X_nonlinear_k, U_ref_k, dt_test, P);
        
        % Track error between methods
        errors_sim(k) = norm(X_linear_sim(:, k+1) - X_nonlinear_sim(:, k+1));
        
        % Additional analysis: Compare linearization quality at each step
        if k <= 5  % Only for first few steps to avoid clutter
            % What would happen if we linearized around current nonlinear state?
            [A_true, B_true, ~, w_true] = get_linearized_dynamics(X_nonlinear_k, U_ref_k, P);
            X_linear_alt = A_true * X_nonlinear_k + B_true * U_ref_k + w_true;
            
            error_ref_lin = norm(X_linear_sim(:, k+1) - X_nonlinear_sim(:, k+1));
            error_true_lin = norm(X_linear_alt - X_nonlinear_sim(:, k+1));
            
            fprintf('    Step %d: ref_linearization_error=%.2e, true_linearization_error=%.2e\n', ...
                k, error_ref_lin, error_true_lin);
        end
    end
    
    % Create comprehensive comparison figure
    if exist('FigureManager', 'class')
        try
            fm = FigureManager();
            fh_val = fm.newFigure('Linearization_Validation_Comparison');
        catch
            fh_val = figure('Name', 'Linearization_Validation_Comparison');
        end
    else
        fh_val = figure('Name', 'Linearization_Validation_Comparison');
    end
    
    fig_val = fh_val;
    
    TL_val = tiledlayout(fig_val, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(TL_val, 'Linear vs Nonlinear Dynamics Comparison');
    
    % State names for plotting
    state_names = {'x (m)', 'y (m)', 'v_x (m/s)', 'v_y (m/s)', '\theta (rad)', '\omega (rad/s)', 'mass (kg)'};
    
    % Plot first 7 states (skip error plot for now)
    for i = 1:7
        nexttile;
        plot(t_sim, X_linear_sim(i, :), 'b-', 'LineWidth', 2, 'DisplayName', 'Linear');
        hold on;
        plot(t_sim, X_nonlinear_sim(i, :), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Nonlinear');
        grid on;
        ylabel(state_names{i});
        if i == 1
            legend('Location', 'best');
        end
        if i >= 6
            xlabel('Time (s)');
        end
    end
    
    % Add error plot in the 8th tile
    nexttile;
    semilogy(t_sim(2:end), errors_sim(1:length(t_sim)-1), 'k-', 'LineWidth', 2);
    grid on;
    ylabel('||Error|| (log scale)');
    xlabel('Time (s)');
    title('Prediction Error');
    
    % Calculate and display key metrics
    max_error_sim = max(errors_sim);
    avg_error_sim = mean(errors_sim);
    final_error_sim = errors_sim(end);
    
    fprintf('   Extended simulation results:\n');
    fprintf('     Maximum error: %.2e\n', max_error_sim);
    fprintf('     Average error: %.2e\n', avg_error_sim);
    fprintf('     Final error: %.2e\n', final_error_sim);
    fprintf('     Error growth rate: %.2e per step\n', final_error_sim / n_sim_steps);
    
    % Detailed diagnostic analysis
    fprintf('\n   DIAGNOSTIC ANALYSIS:\n');
    
    % Check if states remain bounded
    max_linear_state = max(abs(X_linear_sim(:)));
    max_nonlinear_state = max(abs(X_nonlinear_sim(:)));
    fprintf('     Max linear state magnitude: %.2e\n', max_linear_state);
    fprintf('     Max nonlinear state magnitude: %.2e\n', max_nonlinear_state);
    
    % Check for NaN or Inf
    if any(~isfinite(X_linear_sim(:)))
        fprintf('     WARNING: Linear simulation contains NaN/Inf values!\n');
    end
    if any(~isfinite(X_nonlinear_sim(:)))
        fprintf('     WARNING: Nonlinear simulation contains NaN/Inf values!\n');
    end
    
    % Analyze error growth pattern
    if length(errors_sim) > 3
        log_errors = log(max(errors_sim, 1e-16));  % Avoid log(0)
        error_growth_rate = (log_errors(end) - log_errors(1)) / length(log_errors);
        fprintf('     Log error growth rate: %.3f per step\n', error_growth_rate);
        
        if error_growth_rate > 0.1
            fprintf('     ERROR: Exponential error growth detected!\n');
            fprintf('     This suggests:\n');
            fprintf('       - Unstable linearized system\n');
            fprintf('       - Time step too large\n');
            fprintf('       - Linearization errors in get_linearized_dynamics\n');
        end
    end
    
    % Check linearization matrices for stability
    X_check = X_ref(:, 1);
    U_check = U_ref(:, 1);
    [A_check, ~, ~, ~] = get_linearized_dynamics(X_check, U_check, P);
    A_discrete = eye(7) + A_check * dt_test;  % Forward Euler discretization
    eigenvals = eig(A_discrete);
    max_eigenval = max(abs(eigenvals));
    fprintf('     Discrete A matrix max eigenvalue: %.3f\n', max_eigenval);
    
    if max_eigenval > 1.0
        fprintf('     WARNING: Unstable discrete system (max eigenvalue > 1)\n');
        fprintf('     Eigenvalues: ');
        fprintf('%.3f ', abs(eigenvals));
        fprintf('\n');
        
        if max_eigenval > 1.1
            fprintf('     CRITICAL: Highly unstable - reduce time step!\n');
        end
    end
    
    % Additional analysis: compute speeds, angles for comparison  
    speed_linear = sqrt(X_linear_sim(3,:).^2 + X_linear_sim(4,:).^2);
    speed_nonlinear = sqrt(X_nonlinear_sim(3,:).^2 + X_nonlinear_sim(4,:).^2);
    
    % Create second detailed comparison figure
    if exist('FigureManager', 'class')
        try
            fh_val2 = fm.newFigure('Linearization_Detailed_States');
        catch
            fh_val2 = figure('Name', 'Linearization_Detailed_States');
        end
    else
        fh_val2 = figure('Name', 'Linearization_Detailed_States');
    end
    
    fig_val2 = fh_val2;
    
    TL_val2 = tiledlayout(fig_val2, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(TL_val2, 'Key Derived States Comparison (Linear vs Nonlinear)');
    
    % 1. Speed comparison
    nexttile;
    plot(t_sim, speed_linear, 'b-', 'LineWidth', 2, 'DisplayName', 'Linear');
    hold on;
    plot(t_sim, speed_nonlinear, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Nonlinear');
    grid on;
    ylabel('Speed Norm (m/s)');
    legend('Location', 'best');
    
    % 2. Mass comparison
    nexttile;
    plot(t_sim, X_linear_sim(7,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Linear');
    hold on;
    plot(t_sim, X_nonlinear_sim(7,:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Nonlinear');
    grid on;
    ylabel('Mass (kg)');
    
    % 3. Theta comparison (convert to degrees)
    nexttile;
    plot(t_sim, rad2deg(X_linear_sim(5,:)), 'b-', 'LineWidth', 2, 'DisplayName', 'Linear');
    hold on;
    plot(t_sim, rad2deg(X_nonlinear_sim(5,:)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Nonlinear');
    grid on;
    ylabel('\theta (deg)');
    
    % 4. Error evolution
    nexttile;
    semilogy(t_sim(2:end), errors_sim(1:length(t_sim)-1), 'k-', 'LineWidth', 2);
    hold on;
    % Add trend line for exponential growth detection
    if length(errors_sim) > 3
        log_errors = log(errors_sim(1:end));
        p = polyfit(t_sim(2:length(log_errors)+1), log_errors, 1);
        trend_errors = exp(polyval(p, t_sim(2:end)));
        plot(t_sim(2:end), trend_errors, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Trend');
        fprintf('     Error growth exponent: %.3f (negative is good)\n', p(1));
        legend('Actual Error', 'Exponential Trend', 'Location', 'best');
    end
    grid on;
    ylabel('||Error|| (log scale)');
    xlabel('Time (s)');
    
    fprintf('   Figures saved for detailed analysis.\n');
end

%% Test 2: Numerical Derivative Validation
fprintf('\n2. NUMERICAL DERIVATIVE VALIDATION:\n');

% Test at the first reference point
X_test = X_ref(:, 1);
U_test = U_ref(:, 1);

% Get analytical derivatives
[A_analytical, B_analytical, ~, ~] = get_linearized_dynamics(X_test, U_test, P);

% Compute numerical derivatives
fprintf('  Computing numerical derivatives...\n');
A_numerical = zeros(size(A_analytical));
B_numerical = zeros(size(B_analytical));

% Numerical partial derivatives w.r.t. state
for i = 1:length(X_test)
    X_plus = X_test; X_plus(i) = X_plus(i) + eps;
    X_minus = X_test; X_minus(i) = X_minus(i) - eps;
    
    f_plus = compute_dynamics_rhs(X_plus, U_test, P);
    f_minus = compute_dynamics_rhs(X_minus, U_test, P);
    
    A_numerical(:, i) = (f_plus - f_minus) / (2 * eps);
end

% Numerical partial derivatives w.r.t. control
for i = 1:length(U_test)
    U_plus = U_test; U_plus(i) = U_plus(i) + eps;
    U_minus = U_test; U_minus(i) = U_minus(i) - eps;
    
    f_plus = compute_dynamics_rhs(X_test, U_plus, P);
    f_minus = compute_dynamics_rhs(X_test, U_minus, P);
    
    B_numerical(:, i) = (f_plus - f_minus) / (2 * eps);
end

% Compare analytical vs numerical
A_error = norm(A_analytical - A_numerical, 'fro');
B_error = norm(B_analytical - B_numerical, 'fro');

fprintf('  A matrix error (Frobenius norm): %.2e\n', A_error);
fprintf('  B matrix error (Frobenius norm): %.2e\n', B_error);

% Find largest individual errors
[max_A_error, max_A_idx] = max(abs(A_analytical(:) - A_numerical(:)));
[max_B_error, max_B_idx] = max(abs(B_analytical(:) - B_numerical(:)));

fprintf('  Largest A error: %.2e at element (%d,%d)\n', max_A_error, ...
    mod(max_A_idx-1, size(A_analytical,1))+1, ceil(max_A_idx/size(A_analytical,1)));
fprintf('  Largest B error: %.2e at element (%d,%d)\n', max_B_error, ...
    mod(max_B_idx-1, size(B_analytical,1))+1, ceil(max_B_idx/size(B_analytical,1)));

%% Test 3: Coordinate Transformation Verification
fprintf('\n3. COORDINATE TRANSFORMATION VERIFICATION:\n');

% Test coordinate transformations at several angles
test_angles = [0, pi/6, pi/4, pi/3, pi/2];
for theta = test_angles
    % Test rotation matrices
    T_I_B = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    T_B_I = [cos(theta), sin(theta); -sin(theta), cos(theta)];
    
    % Should be orthogonal: T_I_B * T_B_I = I
    orthogonal_error = norm(T_I_B * T_B_I - eye(2), 'fro');
    
    if orthogonal_error > 1e-14
        fprintf('  WARNING: Rotation matrix orthogonality error at theta=%.2f: %.2e\n', ...
            theta, orthogonal_error);
    end
    
    % Test thrust vector transformation
    delta_test = pi/8;  % 22.5 degrees
    T_test = 500e3;     % 500 kN
    
    % Manual calculation
    thrust_B = [0; T_test];
    T_B_thrust = [cos(delta_test), -sin(delta_test); sin(delta_test), cos(delta_test)];
    thrust_gimbaled_B = T_B_thrust * thrust_B;
    thrust_I_manual = T_I_B * thrust_gimbaled_B;
    
    % Using the formula from the code
    thrust_I_code = T_test * [cos(theta)*sin(delta_test) - sin(theta)*cos(delta_test); 
                              sin(theta)*sin(delta_test) + cos(theta)*cos(delta_test)];
    
    transform_error = norm(thrust_I_manual - thrust_I_code);
    if transform_error > 1e-10
        fprintf('  WARNING: Thrust transformation error at theta=%.2f, delta=%.2f: %.2e\n', ...
            theta, delta_test, transform_error);
    end
end

%% Test 4: Step Size Sensitivity
fprintf('\n4. STEP SIZE SENSITIVITY ANALYSIS:\n');

dt_values = [0.01, 0.05, 0.1, 0.2, 0.5];
errors_vs_dt = zeros(size(dt_values));

X_base = X_ref(:, 1);
U_base = U_ref(:, 1);

for i = 1:length(dt_values)
    dt = dt_values(i);
    
    % Get linearized dynamics for this dt
    P_temp = P; P_temp.dt = dt;
    [A_dt, B_dt, S_dt, w_dt] = get_linearized_dynamics(X_base, U_base, P_temp);
    
    % Predict using linearized model
    X_pred_linear = A_dt * X_base + B_dt * U_base + w_dt;
    
    % True nonlinear step
    X_pred_nonlinear = simulate_step_2d(X_base, U_base, dt, P);
    
    errors_vs_dt(i) = norm(X_pred_linear - X_pred_nonlinear);
    fprintf('  dt = %.3f: error = %.2e\n', dt, errors_vs_dt(i));
end

% Check if error grows linearly with dt (expected for good linearization)
dt_ratios = dt_values(2:end) ./ dt_values(1:end-1);
error_ratios = errors_vs_dt(2:end) ./ errors_vs_dt(1:end-1);

fprintf('  Expected error scaling with dt:\n');
for i = 1:length(dt_ratios)
    fprintf('    dt ratio %.2f -> error ratio %.2f (expected ~%.2f)\n', ...
        dt_ratios(i), error_ratios(i), dt_ratios(i));
end

%% Summary and Recommendations
fprintf('\n=== VALIDATION SUMMARY ===\n');
fprintf('Maximum linearization error: %.2e\n', validation_results.max_error);
fprintf('Average linearization error: %.2e\n', validation_results.avg_error);

if validation_results.max_error > 1e-3
    fprintf('*** LINEARIZATION ACCURACY ISSUES DETECTED ***\n');
    fprintf('Recommendations:\n');
    
    if A_error > 1e-6 || B_error > 1e-6
        fprintf('- Check analytical partial derivatives in get_linearized_dynamics\n');
    end
    
    if any(errors_vs_dt > 0.1)
        fprintf('- Reduce time step size for better linearization accuracy\n');
    end
    
    fprintf('- Verify coordinate transformation math\n');
    fprintf('- Check reference trajectory quality\n');
else
    fprintf('Linearization appears accurate.\n');
end

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


function f = compute_dynamics_rhs(X, U, P)
%COMPUTE_DYNAMICS_RHS Computes right-hand side of dynamics: X_dot = f(X,U)
% This should match the continuous-time dynamics before discretization

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