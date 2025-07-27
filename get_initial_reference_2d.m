function [X_ref, U_ref] = get_initial_reference_2d(P, prev_sol, dt_applied)
%GET_INITIAL_REFERENCE_2D Creates a crude initial guess for the trajectory.
%
% The trajectory is initialized with a simple linear interpolation for states
% and a constant control profile. A good initial guess is not required for
% SCP to work, but can speed up convergence.
%
% INPUTS:
%   P           - SCP parameters structure
%   prev_sol    - Previous SCP solution (optional)
%   dt_applied  - Actual timestep applied since previous SCP call (optional)

N = P.N;
n_x = P.n_states;
n_u = P.n_controls;

X_ref = zeros(n_x, N + 1);
U_ref = zeros(n_u, N);

% Check if we have a previous solution to warm start from
if nargin >= 2 && ~isempty(prev_sol)
    [X_ref, U_ref] = warm_start_from_previous(P, prev_sol, dt_applied);
    return;
end

% --- State Trajectory Initialization ---
% Linearly interpolate all states from initial to target
target_state = [P.x_target; P.y_target; P.vx_target; P.vy_target; ...
                P.theta_target; P.omega_target; P.m_dry + 0.1*P.fuel_mass];

for i = 1:n_x
    X_ref(i, :) = linspace(P.X0(i), target_state(i), N + 1);
end
% Ensure mass doesn't go below dry mass
X_ref(7,:) = max(X_ref(7,:), P.m_dry);

% --- Control Trajectory Initialization ---
% Initialize with realistic thrust profile and gimbal for lateral maneuvering
avg_mass = mean(X_ref(7,:));

% Start with higher thrust for deceleration, then reduce
thrust_profile = linspace(1.2, 0.8, N); % Thrust multiplier profile
hover_thrust = avg_mass * P.g0;
for k = 1:N
    T_init = thrust_profile(k) * hover_thrust;
    U_ref(1, k) = min(max(T_init, P.T_min), P.T_max);
end

% Estimate required gimbal angle for lateral maneuvering
% Simple proportional control to reduce horizontal position error
for k = 1:N
    x_error = X_ref(1, k);  % Horizontal position error
    % Simple proportional gimbal (small angle approximation)
    delta_est = -0.1 * x_error / 1000;  % Scale factor
    U_ref(2, k) = max(-P.delta_max, min(P.delta_max, delta_est));
end

% fprintf('Initial reference: T=%.1f-%.1f kN, gimbal=%.2f-%.2f deg\n', ...
%     min(U_ref(1,:))/1e3, max(U_ref(1,:))/1e3, ...
%     rad2deg(min(U_ref(2,:))), rad2deg(max(U_ref(2,:))));

end

function [X_ref, U_ref] = warm_start_from_previous(P, prev_sol, dt_applied)
%WARM_START_FROM_PREVIOUS Shifts previous solution to create warm start
%
% This function takes the previous SCP solution and shifts it forward in time
% to create a good initial guess for the next SCP iteration.
%
% INPUTS:
%   P           - Current SCP parameters structure
%   prev_sol    - Previous SCP solution structure
%   dt_applied  - Actual time elapsed since previous SCP call (s)

N = P.N;
n_x = P.n_states;
n_u = P.n_controls;

X_ref = zeros(n_x, N + 1);
U_ref = zeros(n_u, N);

% Safe trajectory resampling with error handling
try
    % Get previous solution dimensions and timestep
    N_prev = size(prev_sol.X, 2) - 1;  % Number of control intervals in previous solution
    dt_prev = prev_sol.P_scp.dt;
    
    % Validate inputs and detect mode transitions
    if nargin < 3 || isempty(dt_applied) || dt_applied <= 0
        dt_applied = dt_prev; % Fallback to previous timestep if not provided
    end
    
    % Detect timestep mode transitions (normal ↔ fine computation)
    timestep_ratio = dt_prev / P.dt;
    mode_transition = abs(timestep_ratio - 1.0) > 0.1; % Significant timestep change
    
    if mode_transition
        % Special handling for mode transitions to ensure smoother reference
        time_shift = dt_applied * 0.8; % More conservative time shift for mode transitions
    else
        time_shift = dt_applied; % Normal time shift
    end
    
    % Create time grids for resampling
    t_prev = 0:dt_prev:(N_prev*dt_prev);              % Previous time grid
    t_current = 0:P.dt:(N*P.dt);                      % Current time grid
    
    % Shift previous time grid to account for elapsed time
    t_prev_shifted = t_prev + time_shift;
    
    % Validate dimensions and time coverage
    if length(t_prev_shifted) ~= size(prev_sol.X, 2)
        error('Previous solution state grid size mismatch: t_grid=%d, X_size=%d', ...
              length(t_prev_shifted), size(prev_sol.X, 2));
    end
    
    if N_prev > 0 && size(prev_sol.U, 2) ~= N_prev
        error('Previous solution control grid size mismatch: N_prev=%d, U_size=%d', ...
              N_prev, size(prev_sol.U, 2));
    end
    
    % Check time coverage - need some overlap for reliable interpolation
    min_overlap = 0.5 * max(t_current); % Require at least 50% coverage
    if max(t_prev_shifted) < min_overlap
        error('Insufficient time coverage: prev_max=%.3f, required=%.3f', ...
              max(t_prev_shifted), min_overlap);
    end
    
    % Resample states onto new time grid using robust interpolation
    for i = 1:n_x
        if max(t_prev_shifted) >= max(t_current)
            % Full coverage - use linear interpolation
            X_ref(i, :) = interp1(t_prev_shifted, prev_sol.X(i, :), t_current, 'linear', 'extrap');
        else
            % Partial coverage - interpolate available points, extrapolate constant
            idx_valid = t_current <= max(t_prev_shifted);
            X_ref(i, idx_valid) = interp1(t_prev_shifted, prev_sol.X(i, :), t_current(idx_valid), 'linear');
            % Extrapolate with final value for points beyond coverage
            X_ref(i, ~idx_valid) = prev_sol.X(i, end);
        end
    end
    
    % Resample controls with proper handling of discrete nature
    if N_prev > 0
        % Control time grids (controls applied at interval midpoints)
        t_prev_ctrl = t_prev_shifted(1:end-1) + dt_prev/2;
        t_current_ctrl = t_current(1:end-1) + P.dt/2;
        
        for i = 1:n_u
            if max(t_prev_ctrl) >= max(t_current_ctrl)
                % Full coverage - use nearest neighbor for discrete controls
                U_ref(i, :) = interp1(t_prev_ctrl, prev_sol.U(i, :), t_current_ctrl, 'nearest', 'extrap');
            else
                % Partial coverage
                idx_valid = t_current_ctrl <= max(t_prev_ctrl);
                U_ref(i, idx_valid) = interp1(t_prev_ctrl, prev_sol.U(i, :), t_current_ctrl(idx_valid), 'nearest');
                % Extrapolate with final control value
                U_ref(i, ~idx_valid) = prev_sol.U(i, end);
            end
        end
    end
    
    % Ensure first state matches current initial condition exactly
    X_ref(:, 1) = P.X0;
    
    % Validate resampled solution
    if any(isnan(X_ref(:))) || any(isnan(U_ref(:)))
        error('NaN values in resampled trajectory');
    end
    
    % Ensure physical constraints are maintained
    X_ref(7, :) = max(X_ref(7, :), P.m_dry); % Mass constraint
    U_ref(1, :) = max(P.T_min, min(P.T_max, U_ref(1, :))); % Thrust limits
    U_ref(2, :) = max(-P.delta_max, min(P.delta_max, U_ref(2, :))); % Gimbal limits
    
catch ME
    % Fallback to cold start if resampling fails
%     fprintf('  Warning: Warm start failed (%s), using cold start\n', ME.message);
    [X_ref, U_ref] = cold_start_reference(P);
end
end

function [X_ref, U_ref] = cold_start_reference(P)
%COLD_START_REFERENCE Generate initial reference from scratch
N = P.N;
n_x = P.n_states;
n_u = P.n_controls;

X_ref = zeros(n_x, N + 1);
U_ref = zeros(n_u, N);

% State trajectory: linearly interpolate from initial to target
target_state = [P.x_target; P.y_target; P.vx_target; P.vy_target; ...
               P.theta_target; P.omega_target; P.m_dry + 0.1*P.fuel_mass];

for i = 1:n_x
    X_ref(i, :) = linspace(P.X0(i), target_state(i), N + 1);
end
X_ref(7,:) = max(X_ref(7,:), P.m_dry);  % Ensure mass doesn't go below dry mass

% Control trajectory: reasonable initial guess
avg_mass = mean(X_ref(7,:));
thrust_profile = linspace(1.2, 0.8, N);
hover_thrust = avg_mass * P.g0;

for k = 1:N
    T_init = thrust_profile(k) * hover_thrust;
    U_ref(1, k) = min(max(T_init, P.T_min), P.T_max);
    
    % Simple proportional gimbal for lateral maneuvering
    x_error = X_ref(1, k);
    delta_est = -0.1 * x_error / 1000;
    U_ref(2, k) = max(-P.delta_max, min(P.delta_max, delta_est));
end

% fprintf('  Cold start: Generated initial reference\n');
end