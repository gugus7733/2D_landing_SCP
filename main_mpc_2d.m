%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    2D ROCKET SOFT LANDING WITH MPC                      %
%                                                                         %
%   This script simulates the terminal landing phase of a Falcon 9-like   %
%   rocket in 2D using a Model Predictive Control (MPC) strategy based    %
%   on Successive Convexification (SCP).                                  %
%                                                                         %
%   Date: July 20, 2025                                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

fm = FigureManager('exportDir','./figures_2d','visible',true, 'export_figs', 0);

%% Reference Frame & Coordinate System Definitions
%
% --- Inertial Frame (I) ---
% A flat-Earth, non-rotating frame fixed to the landing pad.
% - Origin: (0,0) is the target landing spot.
% - x-axis: Horizontal direction (right is positive).
% - y-axis: Vertical direction (altitude), positive upwards.
% - z-axis: Perpendicular to the x-y plane, pointing out of the page towards the observer.
%
% --- Rocket Body Frame (B) ---
% A frame attached to the rocket's Center of Mass (CoM).
% - xb-axis: Orthogonal to the rocket's longitudinal axis, going right when the rocket is upright.
% - yb-axis: Along the rocket's longitudinal axis, pointing from tail to nose.
% - zb-axis: equal to intertial z-axis, pointing out of the page.
% When the rocket is upright, xb points to the right, and yb is up. In this situation, this reference orientation is equal to the inertial frame, thus theta = 0.
% A Theta of 180° means the rocket is upside down, with xb pointing left and yb poiting down.
%
% --- Thrust Vector (T) ---
% A vector that represents the thrust force applied by the rocket's engines.
% It is directed along the yb-axis, rotated by a gimbal angle delta (δ) along the z-axis.
% A positive gimbal angle gives a thrust vector that points towards the negative inertial x-axis, and gives a negative moment around the zb-axis which make the rocket have a negative omega (=dtheta).
% thrust_vector = [0;thrust_amount]
%
% --- Aero forces vector (T) ---
% Opposed to the speed vector in the reference frame. Defines the angle of attack (α) and is used to compute the aerodynamic forces.
%
%
%
% --- Key Angles ---
% - theta (θ): Rocket pitch angle. The angle of the rocket's yb-axis
%   relative to the inertial y-axis. theta = 0 means the rocket is
%   perfectly vertical, nose up. Positive theta pitches the nose towards
%   the negative inertial x-axis (rotates counter-clockwise).
% - delta (δ): Engine gimbal angle. The angle of the thrust vector
%   relative to the rocket's positive yb-axis.
%   In order to compute the thrust in the inertial frame for integration, we need to take into account the gimbal:
%   We have xy_thrust_vector = [0;thrust_amount]. In order to express this in the inertial frame, we need to rotate it by the angle +-theta and +-delta (rotation matrices):
%   Which is thrust_vector_I = T_I_B * T_B_thrust * thrust_vector.
% - alpha (α): Angle of attack. The angle between the rocket's aerodynamic forces and the yb axis. The aerodynamic forces are opposed to the speed vector.
%   A zero alpha means the rocket is flying straight, with no lateral aerodynamic forces acting on it (only longitudinal). As the aero center is behind the CoM, a positive alpha will create a negative moment around zb, which stabilizes the rocket.
%
% --- Rotation Matrix (Body to Inertial) ---
% This matrix transforms a vector from the Body frame to the Inertial frame. Let's name the rotation matrices T_A_B which transforms a vector from frame B to frame A.
% Thus, V_A = T_A_B * V_B, V_C = T_C_B * T_B_A * V_A, etc.
% T_B_I = [  cos(θ)   sin(θ) ]
%         [ -sin(θ)   cos(θ) ]
% T_I_B = T_B_I';
% A vector v_I in the inertial frame is found by: v_I = T_I_B * v_B. And T_A_B = T_B_A'.
% v_aero = - v_speed / norm(v_speed); then we express it in the body frame: v_aero_body = T_B_I * v_aero; and alpha = atan2(v_aero_body(1), v_aero_body(2));
% The aero forces apply forces and moments on the rocket in the body frame, so we need to express them in the inertial frame for integration.
% For the thrust vector, we have:
% thrust_vector_B = T_B_thrust * [0;T], with T_B_thrust = [cos(delta) -sin(delta); sin(delta) cos(delta)];

%% Rotation Matrix Functions
% Define rotation matrices as functions for consistent use throughout the code

% Rotation matrix from Body frame to Inertial frame
T_I_B = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];

% Rotation matrix from Inertial frame to Body frame  
T_B_I = @(theta) [cos(theta), sin(theta); -sin(theta), cos(theta)];

% Thrust gimbal rotation matrix in Body frame
T_B_thrust = @(delta) [cos(delta), -sin(delta); sin(delta), cos(delta)];


%% Physical & Model Parameters
P = struct();

% === Simulation & Horizon ===
P.dt_sim        = 0.01;              % Simulation timestep (s)
P.N_scp_steps_open_loop = 30;         % Number of dt_scp steps to apply in open loop before recomputing SCP (0 = classic MPC)
P.dt_scp        = 0.2;               % SCP optimization timestep (s)
P.max_iters_scp = 10;
P.fine_computation_time = 0.0;     % s, time before touchdown to switch to fine computation
P.fine_computation_dt = 0.1;      % s, fine SCP dt (much smaller timestep)
P.fine_computation_n_iter = 3;     % fine SCP max iterations (more iterations for accuracy)
P.T_max_mission = 50;               % Max mission time (s)
P.min_horizon_time = 2*P.fine_computation_dt;            % Minimum SCP horizon time (s)

% === Vehicle Properties ===
P.m_dry         = 22000;             % Dry mass (kg)
P.m0            = 35000;             % Initial wet mass (kg)
P.fuel_mass     = P.m0 - P.m_dry;    % Available propellant (kg)
P.L_rocket      = 40;                % Approx. length of rocket (m)
% Moment of inertia (Iyy) modeled as a slender rod, varies with mass
% I = (1/12)*m*L^2. We define a function for it later.
P.Iyy_func      = @(m) (1/12) * m * P.L_rocket^2;
P.L_com_from_base = 0.4 * P.L_rocket; % Center of Mass location from base (m)
P.L_cop_from_base = 0.8 * P.L_rocket; % Center of Pressure from base (m)

% === Initial State [x, y, vx, vy, theta, omega, m]' ===
P.x0      = -1000;              % m, initial horizontal position
P.y0      = 5000;              % m, initial altitude
P.vx0     = 100;              % m/s, initial horizontal velocity
P.vy0     = -200;              % m/s, initial vertical velocity
P.alpha0  = 15*d2r;            % rad, initial AoA
P.theta0  = computeInitialSlope(P.alpha0, [P.vx0; P.vy0]);        % rad, initial pitch
P.omega0  = deg2rad(-20);        % rad/s, initial angular velocity

% === Target State ===
P.x_target     = 0;
P.y_target     = 0;
P.vx_target    = 0;
P.vy_target    = 0;
P.theta_target = 0;
P.omega_target = 0;

% === Propulsion ===
P.T_min       = 0.2 * 934e3;       % N, minimum thrust (approx. 1 Merlin engine)
P.T_max       = 1e6;             % N, maximum thrust
P.delta_max   = deg2rad(30.0);      % rad, max gimbal angle
P.Isp         = 282;               % s, specific impulse

% === Aerodynamics ===
P.Cd          = 0.8;               % Axial drag coefficient (approx.)
P.C_N_alpha   = 5;               % Normal force coeff. derivative (per rad)
P.A_ref       = 10.8;              % m^2, reference area (diameter 3.7m)
P.C_damp      = 1e7;               % Aerodynamic damping coefficient (set >0 to enable damping)

% === Environment ===
P.g0          = 9.80665;           % m/s^2, standard gravity

% === Cost Weights & Penalties ===
P.w_T         = 1e-10;              % Thrust magnitude penalty
P.w_delta     = 0;              % Gimbal angle penalty
P.w_dT        = 0;              % Thrust rate penalty
P.w_ddelta    = 1e-1;               % Gimbal rate penalty
P.w_omega     = 1e-4;               % Angular velocity penalty

% === Trust Region Parameters (SCP Adaptive) ===
% trust regions init
P.trust_init_T     = 1 * P.T_max;   % Initial trust in thrust change
P.trust_init_delta = 1 * P.delta_max; % Initial trust in gimbal change
P.trust_init_vx    = 200;               % m/s trust for vx changes
P.trust_init_vy    = 200;               % m/s trust for vy changes  
P.trust_init_omega = deg2rad(90);      % rad/s trust for omega changes
% trust regions bounds
P.trust_min_T      = 0.5 * P.T_max;   % Minimum trust region for thrust
P.trust_max_T      = 1 * P.T_max;    % Maximum trust region for thrust
P.trust_min_delta  = 0.5 * P.delta_max; % Minimum trust region for gimbal
P.trust_max_delta  = 1.0 * P.delta_max;  % Maximum trust region for gimbal
P.trust_min_vx     = 50;                % minimum vx trust region
P.trust_max_vx     = 200;              % maximum vx trust region
P.trust_min_vy     = 50;                % minimum vy trust region
P.trust_max_vy     = 200;              % maximum vy trust region
P.trust_min_omega  = deg2rad(10);       % minimum omega trust
P.trust_max_omega  = deg2rad(90);      % maximum omega trust

% Trust region adaptation
P.trust_shrink     = 0.9;              % Trust region shrink factor
P.trust_expand     = 1.25;             % Trust region expand factor
P.slack_trigger    = 1e-3;             % Slack threshold for trust region adaptation

% === Enhanced Slack Management Parameters ===
P.w_slack_initial = 1e-6;           % Initial slack weight (lower for early convergence)
P.w_slack_terminal = 1e6;           % Terminal slack weight (high to force zero slack)
P.slack_decay_alpha = 2.0;          % Slack weight progression exponent
P.slack_decay_beta = 3.0;           % Time progression exponent (higher = sharper transition)
P.slack_max_initial = 400;          % Initial maximum allowed slack magnitude (m/s, rad/s)
P.slack_max_terminal = 1e-1;        % Terminal maximum allowed slack magnitude
P.slack_retry_multiplier = 100;     % Factor to increase slack bounds on retry
P.max_slack_retries = 3;            % Maximum number of retry attempts with relaxed slack and trust

% === SCP Tolerances ===
P.tol_cost    = 1e-12;              % Cost change tolerance for convergence
P.tol_slack   = 1e-12;              % Slack tolerance for convergence
P.tol_slack_initial = 1e-6;         % Initial slack tolerance (more relaxed)
P.tol_slack_terminal = 1e-12;       % Terminal slack tolerance (very strict)

% === Legacy parameters (kept for compatibility) ===
P.trust_T     = 0.3 * P.T_max;     % Trust region for thrust
P.trust_delta = deg2rad(10);       % Trust region for gimbal

%% Debug & Validation Options : true or false
isDebug = false;
n_debug_runs = 1;
validate_linearization_flag = false; % Set to true to validate linearization accuracy

%% Control Replay System Parameters
replay_control = false;           % Set to true to enable control replay mode
t_replay_control = 25;            % Time until which replay is active (s)
control_replay_filename = 'control_replay.mat'; % Filename for saved control log

%% Open-Loop Control Verification Parameters
open_loop = false;                % Set to true to enable open-loop verification mode
t_open_loop = 1;              % Time to switch to open-loop control (s)
open_loop_horizon_limit = [];     % Max commands to use (empty = all available)


%% MPC Simulation Loop
% Initialize simulation state
current_state = [P.x0; P.y0; P.vx0; P.vy0; P.theta0; P.omega0; P.m0];
sim_history = struct('t', [], 'X', [], 'U', [], 'vectors', struct('velocity', [], 'orientation', [], 'thrust', []));
last_scp_sol = [];
last_scp_log = [];
scp_debug_history = struct('time', [], 'log', [], 'sol', []);
total_scp_calls = 0;
sim_time = 0;

% Initialize control logging for crash recovery
control_log = control_replay_utils('initialize');

% Initialize open-loop control system
open_loop_active = false;
open_loop_control_sequence = [];

% N-step open-loop execution state variables
n_step_open_loop_active = false;
n_step_control_sequence = [];
n_step_dt_sequence = [];
n_step_command_index = 1;
n_step_remaining_steps = 0;
last_scp_call_time = 0;
open_loop_dt_sequence = [];
open_loop_command_index = 1;
open_loop_predicted_trajectory = [];
open_loop_start_time = [];
open_loop_achieved_trajectory = struct('t', [], 'X', []);

% Initialize replay control system if enabled
replay_log = [];
if replay_control
    replay_log = control_replay_utils('load_replay', control_replay_filename);
    if ~isempty(replay_log)
        [is_valid, max_time_available] = control_replay_utils('validate_replay_coverage', replay_log, t_replay_control);
        if is_valid
            fprintf('Replay mode enabled: using controls until t=%.2f s\n', t_replay_control);
        else
            fprintf('Warning: Replay file only covers %.2f s, but t_replay_control=%.2f s\n', max_time_available, t_replay_control);
            fprintf('Replay will be used until %.2f s, then switch to SCP\n', max_time_available);
            t_replay_control = max_time_available;
        end
    else
        fprintf('Warning: Replay file %s not found. Disabling replay mode.\n', control_replay_filename);
        replay_control = false;
    end
end

% Validate enhanced slack management parameters
[is_valid, error_msg] = slack_management_utils('validate_parameters', P);
if ~is_valid
    warning('Slack management parameter validation failed: %s. Using legacy slack system.', error_msg);
end

fprintf('=== 2D MPC Simulation Started ===\n');
fprintf('Initial State: x=%.1f, y=%.1f, vx=%.1f, vy=%.1f, th=%.1f deg\n', ...
    current_state(1), current_state(2), current_state(3), current_state(4), rad2deg(current_state(5)));

if is_valid
    fprintf('Enhanced slack management: w_slack %.1e -> %.1e, s_max %.1f -> %.1e\n', ...
        P.w_slack_initial, P.w_slack_terminal, P.slack_max_initial, P.slack_max_terminal);
end

try
    while sim_time < P.T_max_mission
        % --- Landing Check ---
        if current_state(2) < 0.2 % Altitude check (y)
            fprintf('\n>>> Landing detected at t=%.2f s <<<\n', sim_time);
            break;
        end

        % --- Determine Control: Replay, Open-Loop Verification, N-Step Open-Loop, or SCP ---
        if replay_control && sim_time < t_replay_control
            % Use replay control
            [T_applied, delta_applied] = control_replay_utils('get_replay_control', replay_log, sim_time);
            fprintf('\nt=%.2f s: Using replay control (T=%.0f N, delta=%.2f deg)\n', ...
                sim_time, T_applied, rad2deg(delta_applied));
            dt_scp_current = P.dt_scp; % Use standard timestep for replay
        elseif open_loop_active
            % Use open-loop verification control from captured sequence
            if open_loop_command_index <= size(open_loop_control_sequence, 2)
                T_applied = open_loop_control_sequence(1, open_loop_command_index);
                delta_applied = open_loop_control_sequence(2, open_loop_command_index);
                dt_scp_current = open_loop_dt_sequence(open_loop_command_index);
                fprintf('\nt=%.2f s: Open-loop verification [%d/%d] (T=%.0f N, delta=%.2f deg, dt=%.3f s)\n', ...
                    sim_time, open_loop_command_index, size(open_loop_control_sequence, 2), ...
                    T_applied, rad2deg(delta_applied), dt_scp_current);
                open_loop_command_index = open_loop_command_index + 1;
            else
                % No more open-loop commands available
                fprintf('\n*** OPEN-LOOP COMMANDS EXHAUSTED at t=%.2f s ***\n', sim_time);
                fprintf('Applied %d commands from captured sequence\n', size(open_loop_control_sequence, 2));
                break; % End simulation
            end
        elseif n_step_open_loop_active && n_step_remaining_steps > 0
            % Use N-step open-loop control from captured SCP sequence
            T_applied = n_step_control_sequence(1, n_step_command_index);
            delta_applied = n_step_control_sequence(2, n_step_command_index);
            dt_scp_current = n_step_dt_sequence(n_step_command_index);
            
            % Update N-step counters
            current_step = n_step_command_index - 1; % For logging (1-indexed)
            total_steps = size(n_step_control_sequence, 2);
            n_step_command_index = n_step_command_index + 1;
            n_step_remaining_steps = n_step_remaining_steps - 1;
            
            % Check if N-step sequence is exhausted
            if n_step_remaining_steps == 0
                n_step_open_loop_active = false;
            end
        else
            % Use SCP optimization
            if replay_control && sim_time >= t_replay_control
                fprintf('\n*** REPLAY MODE ENDED - SWITCHING TO SCP ***\n');
                replay_control = false; % Disable replay for remaining simulation
            end
            
            % Check if we should switch to open-loop mode
            if open_loop && sim_time >= t_open_loop && ~open_loop_active
                fprintf('\n*** SWITCHING TO OPEN-LOOP MODE at t=%.2f s ***\n', sim_time);
                if ~isempty(last_scp_sol)
                    % Capture the current SCP solution for open-loop execution
                    open_loop_control_sequence = last_scp_sol.U;
                    
                    % Store the predicted trajectory for comparison
                    open_loop_predicted_trajectory = last_scp_sol.X;
                    open_loop_start_time = sim_time;
                    
                    % Determine timesteps for each command based on current SCP mode
                    if T_est <= P.fine_computation_time
                        dt_for_sequence = P.fine_computation_dt;
                    else
                        dt_for_sequence = P.dt_scp;
                    end
                    
                    % Create timestep sequence
                    n_commands = size(open_loop_control_sequence, 2);
                    if isempty(open_loop_horizon_limit)
                        commands_to_use = n_commands;
                    else
                        commands_to_use = min(n_commands, open_loop_horizon_limit);
                    end
                    
                    open_loop_control_sequence = open_loop_control_sequence(:, 1:commands_to_use);
                    open_loop_dt_sequence = repmat(dt_for_sequence, 1, commands_to_use);
                    open_loop_predicted_trajectory = open_loop_predicted_trajectory(:, 1:commands_to_use+1);
                    
                    fprintf('Captured SCP solution: %d commands, dt=%.3f s each\n', ...
                            commands_to_use, dt_for_sequence);
                    fprintf('Total open-loop horizon: %.2f s\n', commands_to_use * dt_for_sequence);
                    
                    open_loop_active = true;
                    open_loop_command_index = 1;
                    
                    % Initialize achieved trajectory with current state
                    open_loop_achieved_trajectory.t = sim_time;
                    open_loop_achieved_trajectory.X = current_state;
                    
                    % Apply first open-loop command immediately
                    T_applied = open_loop_control_sequence(1, 1);
                    delta_applied = open_loop_control_sequence(2, 1);
                    dt_scp_current = open_loop_dt_sequence(1);
                    fprintf('t=%.2f s: Open-loop control [1/%d] (T=%.0f N, delta=%.2f deg, dt=%.3f s)\n', ...
                        sim_time, commands_to_use, T_applied, rad2deg(delta_applied), dt_scp_current);
                    open_loop_command_index = 2;
                else
                    fprintf('Warning: No SCP solution available for open-loop mode. Continuing with SCP.\n');
                    open_loop = false; % Disable open-loop mode
                end
            end
            
            % Skip SCP if we just switched to open-loop mode
            if open_loop_active && sim_time >= t_open_loop
                % Controls already set above, skip SCP computation
            else
                % Estimate time to touchdown to set optimization horizon
                [T_est, T_horizon] = estimate_time_to_touchdown(current_state, P);
            
                % Detect mode transitions for better reference trajectory handling
                previous_fine_mode = false;
                if ~isempty(last_scp_sol) && isfield(last_scp_sol, 'P_scp')
                    previous_fine_mode = (last_scp_sol.P_scp.dt <= P.fine_computation_dt * 1.1);
                end
                
                % Check if we're in fine computation phase (close to landing)
                if T_est <= P.fine_computation_time
                    % Fine computation: smaller timestep, more iterations
                    dt_scp_current = P.fine_computation_dt;
                    max_iters_current = P.fine_computation_n_iter;
                    
                    if ~previous_fine_mode
                        fprintf('\n*** FINE COMPUTATION MODE ACTIVATED ***\n');
                    end
                    fprintf('t=%.2f s: Fine SCP (dt=%.3f s, max_iters=%d, Horizon=%.2f s)\n', ...
                        sim_time, dt_scp_current, max_iters_current, T_horizon);
                else
                    % Normal computation: standard timestep, fewer iterations  
                    dt_scp_current = P.dt_scp;
                    max_iters_current = P.max_iters_scp;
                    
                    if previous_fine_mode
                        fprintf('\n*** RETURNING TO NORMAL COMPUTATION MODE ***\n');
                    end
                    fprintf('\nt=%.2f s: Normal SCP (dt=%.3f s, max_iters=%d, Horizon=%.2f s)\n', ...
                        sim_time, dt_scp_current, max_iters_current, T_horizon);
                end
                
                N_scp = round(T_horizon / dt_scp_current);
                N_scp = max(N_scp, 3); % Ensure a minimum number of steps
    
                % Run the SCP optimization with adaptive parameters and warm start
                % Pass simulation time for enhanced slack management
    %             if (sim_time == 0)
    %                 [scp_sol, scp_log] = run_scp_2d(current_state, T_horizon, N_scp, P, dt_scp_current, max_iters_current, last_scp_sol, sim_time);
    %             else
    %                 [scp_sol, scp_log] = run_scp_2d(current_state, T_horizon, N_scp, P, dt_scp_current, max_iters_current, last_scp_sol, sim_time);
    %             end
                
                [scp_sol, scp_log] = run_scp_2d(current_state, T_horizon, N_scp, P, dt_scp_current, max_iters_current, last_scp_sol, sim_time);
    
                total_scp_calls = total_scp_calls + 1;
    
                if isempty(scp_sol)
                    scp_sol = last_scp_sol;
                    scp_log = last_scp_log;
                else
                    fprintf('  SCP converged in %d iterations, max retries = %d. Cost: %.2e, Slack: %.2e\n', ...
                        length(scp_log.cost), max(scp_log.retry_level), scp_log.cost(end), scp_log.slack(end));
                end
                
                last_scp_sol = scp_sol;
                last_scp_log = scp_log;
                
                % Validate linearization accuracy (only for first SCP call)
                if validate_linearization_flag && total_scp_calls == 1
                    fprintf('\n--- LINEARIZATION VALIDATION ---\n');
                    validate_linearization(scp_sol.X, scp_sol.U, current_state, T_horizon, N_scp, P);
                    debug_specific_linearization_issues(scp_sol.X, scp_sol.U, current_state, T_horizon, N_scp, P);
                end
    
                    % Apply first optimal control from the SCP solution
                    T_applied = last_scp_sol.U(1,1);
                    delta_applied = last_scp_sol.U(2,1);
                    
                    % Check if should activate N-step open-loop mode
                    % Only activate if not in verification open-loop mode
                    if P.N_scp_steps_open_loop > 0 && ~n_step_open_loop_active && ~open_loop_active
                        % Determine how many steps to capture from SCP solution
                        available_steps = size(last_scp_sol.U, 2);
                        n_steps_to_capture = min(P.N_scp_steps_open_loop, available_steps);
                        
                        % Only activate if more than 1 step available (first step already applied above)
                        if n_steps_to_capture > 1
                            % Capture control sequence
                            n_step_control_sequence = last_scp_sol.U(:, 1:n_steps_to_capture);
                            
                            % Create timestep sequence (timestep-agnostic)
                            % Use consistent timestep from current SCP mode
                            n_step_dt_sequence = repmat(dt_scp_current, 1, n_steps_to_capture);
                            
                            % Initialize N-step execution state
                            n_step_remaining_steps = n_steps_to_capture - 1; % -1 because first step applied immediately
                            n_step_command_index = 2; % Next step to apply (1 is already applied)
                            n_step_open_loop_active = true;
                            last_scp_call_time = sim_time;
                        end
                    end
            end
        end
        
        % Log control for crash recovery
        control_log = control_replay_utils('log_control', control_log, sim_time, T_applied, delta_applied);

        % Simulate forward using the nonlinear model
        % Integrate over one SCP interval (dt_scp_current) using smaller sim steps (dt_sim)
        temp_state = current_state;
        num_sim_steps = round(dt_scp_current / P.dt_sim);

        for i = 1:num_sim_steps
            % Extract state components for vector calculations
            x_pos     = temp_state(1);
            y_pos     = temp_state(2);
            vx        = temp_state(3);
            vy        = temp_state(4);
            theta     = temp_state(5);
            omega     = temp_state(6);
            
            % Calculate vectors for visualization
            % Velocity vector (blue) - inertial frame
            velocity_vector = [vx; vy];
            
            % Rocket orientation vector (black) - body y-axis in inertial frame
            rocket_up_body = [0; 1]; % Unit vector along rocket body y-axis
            orientation_vector = T_I_B(theta) * rocket_up_body;
            
            % Thrust vector (green) - thrust direction in inertial frame
            thrust_body = [0; T_applied]; % Thrust along body y-axis
            thrust_gimbal_body = T_B_thrust(delta_applied) * thrust_body;
            thrust_vector = T_I_B(theta) * thrust_gimbal_body;
            
            % Log state and vectors at each small sim step
            if isempty(sim_history.t)
                sim_history.t = sim_time;
                sim_history.X = temp_state;
                sim_history.U = [T_applied; delta_applied];
                sim_history.vectors.velocity = velocity_vector;
                sim_history.vectors.orientation = orientation_vector;
                sim_history.vectors.thrust = thrust_vector;
            else
                sim_history.t(end+1) = sim_time;
                sim_history.X(:,end+1) = temp_state;
                sim_history.U(:,end+1) = [T_applied; delta_applied];
                sim_history.vectors.velocity(:,end+1) = velocity_vector;
                sim_history.vectors.orientation(:,end+1) = orientation_vector;
                sim_history.vectors.thrust(:,end+1) = thrust_vector;
            end

            % Simulate one small step
            temp_state = simulate_step_2d(temp_state, [T_applied; delta_applied], P.dt_sim, P);
            sim_time = sim_time + P.dt_sim;
        end
        current_state = temp_state;
        
        % Record achieved trajectory during open-loop mode
        if open_loop_active
            open_loop_achieved_trajectory.t(end+1) = sim_time;
            open_loop_achieved_trajectory.X(:,end+1) = current_state;
        end

%         fprintf('  State update: x=%.1f, y=%.1f, vx=%.1f, vy=%.1f, th=%.1f deg\n',...
%             current_state(1), current_state(2), current_state(3), current_state(4), rad2deg(current_state(5)));

        if isDebug
            % Store SCP debug info
            scp_debug_history.time(end+1) = sim_time;
            scp_debug_history.log{end+1} = scp_log;
            scp_debug_history.sol{end+1} = scp_sol;
            if (numel(scp_debug_history.time) >= n_debug_runs)
                break
            end
        end
    end

catch ME
    % Simulation crashed - save trajectory for replay and display error info
    fprintf('\n!!! SIMULATION CRASHED at t=%.2f s !!!\n', sim_time);
    fprintf('Error: %s\n', ME.message);
    fprintf('Saving trajectory to %s for replay...\n', control_replay_filename);
    
    % Save control log for crash replay
    control_replay_utils('save_log', control_log, control_replay_filename);
    
    fprintf('Crash trajectory saved. Set replay_control=true and t_replay_control=%.2f to replay until crash.\n', sim_time);
    
    % Re-throw error to allow debugging if needed
    rethrow(ME);
end

fprintf('\n=== MPC Simulation Complete ===\n');
fprintf('Total time: %.2f s\n', sim_time);
fprintf('Total SCP calls: %d\n', total_scp_calls);

% Open-loop mode summary
if open_loop_active || (open_loop && sim_time >= t_open_loop)
    if ~isempty(open_loop_control_sequence)
        commands_applied = open_loop_command_index - 1;
        total_commands = size(open_loop_control_sequence, 2);
        fprintf('=== OPEN-LOOP VERIFICATION SUMMARY ===\n');
        fprintf('Open-loop started at: t=%.2f s\n', t_open_loop);
        fprintf('Commands applied: %d/%d\n', commands_applied, total_commands);
        if commands_applied < total_commands
            fprintf('Reason for termination: Landing or mission time limit\n');
        else
            fprintf('Reason for termination: All open-loop commands exhausted\n');
        end
        fprintf('Open-loop duration: %.2f s (%.2f s planned)\n', ...
                sim_time - t_open_loop, total_commands * open_loop_dt_sequence(1));
    end
end

fprintf('Final state: x=%.2f, y=%.2f, vx=%.2f, vy=%.2f, th=%.2f, w=%.2f\n',...
    current_state(1), current_state(2), current_state(3), current_state(4),...
    rad2deg(current_state(5)), rad2deg(current_state(6)));
fuel_used = P.m0 - current_state(7);
fprintf('Fuel used: %.1f kg (%.1f%% of available)\n', fuel_used, 100*fuel_used/P.fuel_mass);

% Save control log for future replay
if ~isempty(control_log.t)
    control_replay_utils('save_log', control_log, control_replay_filename);
    fprintf('Control trajectory saved to %s for replay\n', control_replay_filename);
end

%% Visualizations

if (~isDebug)
    visualize_results_2d(sim_history, P, fm);
end

% Interactive trajectory plot
plot_interactive_trajectory_2d(sim_history, P, fm);

% Open-loop trajectory comparison plot
if open_loop_active || (open_loop && sim_time >= t_open_loop)
    if ~isempty(open_loop_predicted_trajectory) && ~isempty(open_loop_achieved_trajectory.X)
        plot_open_loop_comparison(open_loop_predicted_trajectory, open_loop_achieved_trajectory, ...
                                 open_loop_start_time, open_loop_dt_sequence, P, fm);
    end
end

% Add SCP convergence analysis plots
% fprintf('\nGenerating SCP convergence analysis...\n');
for i = 1:length(scp_debug_history.time)
    scp_time = scp_debug_history.time(i);
    scp_log = scp_debug_history.log{i};
    scp_sol = scp_debug_history.sol{i};
    
    % Create convergence plot for each SCP call
    plot_scp_convergence(scp_log, scp_sol, scp_time, fm, P);
    
    % Rename the figures to include SCP call number
    if isfield(fm.handles, 'SCP_Convergence_Analysis')
        fm.handles.(['SCP_Conv_' num2str(i)]) = fm.handles.SCP_Convergence_Analysis;
        fm.handles = rmfield(fm.handles, 'SCP_Convergence_Analysis');
    end
    if isfield(fm.handles, 'SCP_Solution_Analysis')
        fm.handles.(['SCP_Sol_' num2str(i)]) = fm.handles.SCP_Solution_Analysis;
        fm.handles = rmfield(fm.handles, 'SCP_Solution_Analysis');
    end
end

% Create summary plot of all SCP calls
if ~isempty(scp_debug_history.time)
    create_scp_summary_plot(scp_debug_history, fm);
end

fm.exportAll();

%% ====================== Helper Function ==========================================
function [T_est, T_est_relax] = estimate_time_to_touchdown(state, P)
%ESTIMATE_TIME_TO_TOUCHDOWN
%   Estimates realistic and relaxed (padded) time-to-touchdown for 1D descent.
%   Outputs:
%     T_est:      physically realistic estimate (sec)
%     T_est_relax: planning/robustness estimate (sec, always >= T_est)

    y  = state(2);   % altitude (m)
    vy = state(4);   % vertical speed (m/s), negative = down

    margin = 2;

    % --- Physically realistic estimate ---
    if y <= 0
        T_est = 0;
        T_est_relax = 0;
    elseif vy >= -0.5
        % Going up or hovering—use ballistic fallback + margin
        T_est = 2 * sqrt(2 * y / P.g0) + 5;
        T_est_relax = T_est;
    else
        T_est = y / abs(vy);
%         T_est = T_est * 1.25 + 2;  % buffer for braking
        T_est_relax = 1.8*T_est;
    end

    % Enforce limits
    if isfield(P, 'min_horizon_time'), T_est = max(T_est, P.min_horizon_time); end
    if isfield(P, 'T_max_mission'),   T_est = min(T_est, P.T_max_mission);   end

    % Clip relaxed time to global limits
    if isfield(P, 'T_max_mission'),   T_est_relax = min(T_est_relax, P.T_max_mission);   end
end
