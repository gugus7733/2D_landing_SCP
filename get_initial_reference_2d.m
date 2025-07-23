function [X_ref, U_ref] = get_initial_reference_2d(P, prev_sol, dt_applied)
%GET_INITIAL_REFERENCE_2D Creates a crude initial guess for the trajectory.
%
% The trajectory is initialized with a simple linear interpolation for states
% and a constant control profile. A good initial guess is not required for
% SCP to work, but can speed up convergence.

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
    
    % Create time grids for resampling
    t_prev = 0:dt_prev:(N_prev*dt_prev);              % Previous time grid
    t_current = 0:P.dt:(N*P.dt);                      % Current time grid
    
    % Time shift: how much time has elapsed since last SCP call
    time_shift = dt_prev;
    
    % Shift previous time grid to account for elapsed time
    t_prev_shifted = t_prev + time_shift;
    
    % Check if we have enough time coverage for resampling
    if length(t_prev_shifted) > 2 && max(t_prev_shifted) >= max(t_current)
        % Resample states onto new time grid
        for i = 1:n_x
            if length(t_prev_shifted) == size(prev_sol.X, 2)
                X_ref(i, :) = interp1(t_prev_shifted, prev_sol.X(i, :), t_current, 'linear', 'extrap');
            else
                error('Time grid size mismatch');
            end
        end
        
        % Resample controls onto new time grid (use 'nearest' for discrete controls)
        if N_prev > 0
            t_prev_ctrl = t_prev_shifted(1:end-1) + dt_prev/2;  % Control applied at midpoint
            t_current_ctrl = t_current(1:end-1) + P.dt/2;
            
            for i = 1:n_u
                if length(t_prev_ctrl) == size(prev_sol.U, 2)
                    U_ref(i, :) = interp1(t_prev_ctrl, prev_sol.U(i, :), t_current_ctrl, 'nearest', 'extrap');
                else
                    error('Control time grid size mismatch');
                end
            end
        end
        
        % Ensure first state matches current initial condition
        X_ref(:, 1) = P.X0;
        
        % Validate resampled solution
        if any(isnan(X_ref(:))) || any(isnan(U_ref(:)))
            error('NaN values in resampled trajectory');
        end
        
%         fprintf('  Warm start: Resampled previous solution (dt_prev=%.3f → dt_current=%.3f)\n', dt_prev, P.dt);
        
    else
        error('Insufficient time coverage in previous solution');
    end
    
catch ME
    % Fallback to cold start if resampling fails
    fprintf('  Warning: Warm start failed (%s), using cold start\n', ME.message);
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

fprintf('  Cold start: Generated initial reference\n');
end