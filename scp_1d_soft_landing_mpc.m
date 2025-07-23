% simplified 1-D Falcon 9 first-stage terminal landing burn with variable mass.
% -------------------------------------------------------------------------
% This implements a receding horizon (MPC) strategy where:
%   - At each simulation timestep, we run SCP to get optimal trajectory
%   - Apply only the first control action from the SCP solution
%   - Simulate forward using true nonlinear dynamics
%   - Repeat until landing or time runs out
%
% Key features:
%   - Dual timestep: fast simulation dt_sim (0.01s), slower SCP dt_scp (0.5s)
%   - Adaptive horizon: T_final = T_touchdown - T_current
%   - Skip SCP recomputation if states are close to prediction
%   - Reduced verbosity for real-time operation
%
% Author: (Your Name) 2025

%% (0) Housekeeping & Figure Manager -------------------------------------------------
clear; close all; clc

if ~exist('FigureManager','class')
    error('FigureManager class not on path. Add FigureManager before running.');
end
fm = FigureManager('exportDir','./figures','visible',true, 'export_figs', 1);

%% (1) Physical & Model Parameters ---------------------------------------------------

P = struct();

% === Most Used: Timesteps, Windows, Bounds ===
P.dt_sim      = 0.001;             % simulation timestep (s) - small for accurate integration
P.dt_scp      = 0.3;               % fixed SCP optimization timestep (s)
P.T_max_mission = 120;             % maximum mission time (s)
P.min_horizon_time = 2*P.dt_scp;   % minimum horizon time (s)
P.T_min       = 0.2 * 845e3;       % N, minimum thrust
P.T_max       = 845e3;             % N, maximum thrust
P.m_dry       = 15e3;              % dry mass (kg)
P.m0          = 30e3;              % initial mass (kg) (dry + propellant)
P.fuel_mass   = P.m0 - P.m_dry;    % available propellant (kg)
P.p0          = 10000;             % m (altitude above landing pad)
P.v0          = -600;              % m/s (downward negative)
P.p_target    = 0;                 % m
P.v_target    = 0;                 % m/s

% === Fine Computation (Final Phase) ===
P.fine_computation_time = 5.0;     % [s] time before touchdown to switch to fine computation
P.fine_computation_dt = 0.05;      % [s] fine SCP dt
P.fine_computation_n_iter = 100;   % fine SCP max iterations
P.max_iters_fast    = 1;           % fast SCP iterations

% === Model & Environment ===
P.g0          = 9.80665;           % standard gravity (m/s^2)
P.use_var_g   = false;             % toggle for g variation with altitude
P.Re          = 6371000;           % Earth radius (m)
P.mu          = 3.986004418e14;    % Earth GM (m^3/s^2)
P.Isp         = 282;               % s, specific impulse
P.Cd          = 1.5;               % coefficient of drag
P.A_ref       = 10.0;              % m^2 reference area
P.rho0        = 1.225;             % kg/m^3 sea-level density
P.H_scale     = 8500;              % m (exponential scale height)

% === MPC & SCP Settings ===
P.allow_skip_scp = false;          % allow skipping SCP if states are close
P.state_tol_p = 5;                 % position tolerance for skipping SCP recompute (m)
P.state_tol_v = 2;                 % velocity tolerance for skipping SCP recompute (m/s)
P.state_tol_m = 20;                % mass tolerance for skipping SCP recompute (kg)

% === Trust Region Parameters ===
P.trust_init_T = 0.30 * P.T_max;   % large initial trust in thrust change
P.trust_init_v = 150;              % m/s trust for velocity
P.trust_shrink = 0.5;
P.trust_expand = 1.25;

% === Cost Weights & Penalties ===
P.w_T         = 1e-6;              % (placeholder if free-final-time later)
P.w_u         = 1e-6;              % small weight on thrust magnitude^2 (force-based cost)
P.w_s         = 1e6;               % heavy weight on velocity slack (L1)
P.w_du        = 1e-2;              % thrust rate smoothing weight
P.w_term_p    = 0;                 % (hard constraint used)
P.w_term_v    = 0;                 % (hard constraint used)
P.slack_trigger = 1e-3;            % adapt radii threshold

% === Tolerances & Safeguards (Least Used) ===
P.tol_cost    = 1e-6;              % relaxed for MPC
P.tol_slack   = 1e-6;              % relaxed for MPC
P.min_mass    = P.m_dry;           % cannot go below dry mass
P.eps_mass    = 1e-6;
P.eps_rho     = 1e-6;

%% (2) MPC Simulation Loop -----------------------------------------------------------

% Initialize simulation state
sim_time = 0;
current_state = [P.p0; P.v0; P.m0];  % [position; velocity; mass]
sim_history = struct('t', [], 'p', [], 'v', [], 'm', [], 'T_applied', [], 'scp_called', [], 'scp_iters', []);

% SCP solution storage
last_scp_sol = [];
last_scp_time = -inf;
last_predicted_states = [];

% Step counters
scp_step = 0;  % which step in the current SCP solution we're at
total_scp_calls = 0;

% Adaptive iteration variables
current_max_iters = P.max_iters_fast;    % current adaptive max iterations
last_applied_thrust = 0;                % previous applied thrust
low_consistency_count = 0;              % consecutive low consistency count

fprintf('=== MPC Simulation Started ===\n');
fprintf('Simulation dt: %.3f s, SCP dt: %.3f s\n', P.dt_sim, P.dt_scp);
fprintf('Initial state: p=%.1f m, v=%.1f m/s, m=%.1f kg\n', current_state(1), current_state(2), current_state(3));

first_run = true;

while sim_time < P.T_max_mission
    % Check landing condition
    if current_state(1) <= 0  % altitude <= 0
        fprintf('Landing detected at t=%.2f s\n', sim_time);
        break;
    end
    

    % Determine if we need to run SCP
    need_scp = false;
    time_since_scp = sim_time - last_scp_time;

    % Only allow SCP call if at least current_dt_scp has passed since last call, or if no solution yet
    if isempty(last_scp_sol) || isempty(last_predicted_states) || (time_since_scp >= P.dt_scp)
        need_scp = true;
    end

    if ~isempty(last_scp_sol) && ~isempty(last_predicted_states)
        % Check if current state is close to predicted state
        predicted_idx = min(size(last_predicted_states, 1), max(1, round(time_since_scp / P.dt_scp) + 1));
        if predicted_idx <= size(last_predicted_states, 1)
            predicted_state = last_predicted_states(predicted_idx, :)';
            state_error = abs(current_state - predicted_state);
            if state_error(1) < P.state_tol_p && state_error(2) < P.state_tol_v && state_error(3) < P.state_tol_m && P.allow_skip_scp == true
                need_scp = false;
                scp_step = scp_step + 1;
            end
        end
        % Also check if we've used up the current SCP solution
        if scp_step >= length(last_scp_sol.T)
            need_scp = true;
        end
    end
    
    % Run SCP if needed
    if need_scp
        % Calculate remaining time to touchdown (estimate)
        T_remaining = estimateTimeToTouchdown(current_state, P);
        T_remaining = max(T_remaining, P.min_horizon_time);  % ensure minimum horizon

        % Fine computation mode if close to touchdown
        if (first_run == true)
            first_run = false;
            dt_scp = P.dt_scp;
            max_iters = P.fine_computation_n_iter;
        elseif (T_remaining < P.fine_computation_time)
            dt_scp = P.fine_computation_dt;
            max_iters = P.fine_computation_n_iter;
        else
            dt_scp = P.dt_scp;
            max_iters = P.max_iters_fast;
        end

        % Calculate number of SCP intervals using selected timestep
        N_scp = round(T_remaining / dt_scp);
        N_scp = max(N_scp, 2); % at least 2 intervals

        fprintf('t=%.2f: Running SCP (T_rem=%.3f s, N=%d, dt_scp=%.3f, max_iters=%d)\n', ...
            sim_time, T_remaining, N_scp, dt_scp, max_iters);

        % Run SCP optimization with selected iterations
        [scp_sol, scp_log] = runSCP(current_state, T_remaining, N_scp, P, max_iters, dt_scp);
        
        if ~isempty(scp_sol)

            last_scp_sol = scp_sol;
            last_scp_time = sim_time;
            last_predicted_states = [scp_sol.p, scp_sol.v, scp_sol.m];
            scp_step = 1;
            total_scp_calls = total_scp_calls + 1;

            fprintf('  SCP converged in %d iters\n', length(scp_log.cost));
        else
            fprintf('  SCP failed! Using max thrust.\n');
        end
    end
    
    % Apply thrust command
    if ~isempty(last_scp_sol) && scp_step <= length(last_scp_sol.T)
        T_applied = last_scp_sol.T(scp_step);
    else
        % Fallback: emergency thrust
        T_applied = P.T_max;
    end
    
    % Store for next iteration's thrust derivative calculation
    last_applied_thrust = T_applied;
    
    % Simulate one step forward with true nonlinear dynamics
    new_state = simulateOneStep(current_state, T_applied, P.dt_sim, P);
    
    % Log history
    sim_history.t(end+1) = sim_time;
    sim_history.p(end+1) = current_state(1);
    sim_history.v(end+1) = current_state(2);
    sim_history.m(end+1) = current_state(3);
    sim_history.T_applied(end+1) = T_applied;
    sim_history.scp_called(end+1) = need_scp;
    sim_history.scp_iters(end+1) = current_max_iters;
    
    % Update for next iteration
    current_state = new_state;
    sim_time = sim_time + P.dt_sim;
end

fprintf('=== MPC Simulation Complete ===\n');
fprintf('Total time: %.2f s\n', sim_time);
fprintf('Total SCP calls: %d\n', total_scp_calls);
fprintf('Final state: p=%.3f m, v=%.3f m/s, m=%.1f kg\n', ...
    current_state(1), current_state(2), current_state(3));

fuel_used = P.m0 - current_state(3);
fprintf('Fuel used: %.1f kg (%.1f%% of available)\n', fuel_used, 100*fuel_used/P.fuel_mass);

%% (3) Visualization -----------------------------------------------------------------
t_sim = sim_history.t;

% Main trajectory plots
fm.newFigure('MPC_Altitude'); 
plot(t_sim, sim_history.p, 'b-', 'LineWidth', 1.5); 
grid on; xlabel('t (s)'); ylabel('Altitude (m)'); 
title('MPC Trajectory - Altitude');

fm.newFigure('MPC_Velocity'); 
plot(t_sim, sim_history.v, 'r-', 'LineWidth', 1.5); 
grid on; xlabel('t (s)'); ylabel('Velocity (m/s)'); 
title('MPC Trajectory - Velocity');

fm.newFigure('MPC_Mass'); 
plot(t_sim, sim_history.m, 'g-', 'LineWidth', 1.5); 
hold on; yline(P.m_dry, 'k--', 'Dry Mass');
grid on; xlabel('t (s)'); ylabel('Mass (kg)'); 
title('MPC Trajectory - Mass');

fm.newFigure('MPC_Thrust'); 
plot(t_sim, sim_history.T_applied/1e3, 'k-', 'LineWidth', 1.5); 
grid on; xlabel('t (s)'); ylabel('Thrust (kN)'); 
title('MPC Trajectory - Applied Thrust');

% SCP call indicators
fm.newFigure('MPC_SCPCalls');
scp_times = t_sim(sim_history.scp_called == 1);
if ~isempty(scp_times)
    stem(scp_times, ones(size(scp_times)), 'ro', 'LineWidth', 2);
end
grid on; xlabel('t (s)'); ylabel('SCP Called'); 
title('SCP Optimization Calls');
ylim([0, 1.5]);

% Multi-panel summary
afh = fm.newFigure('MPC_Summary'); 
fig = afh.Parent; 
TL = tiledlayout(fig, 5, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; plot(t_sim, sim_history.p, 'LineWidth', 1.3); 
ylabel('Altitude (m)'); grid on; title('MPC Simulation Results');

nexttile; plot(t_sim, sim_history.v, 'LineWidth', 1.3); 
ylabel('Velocity (m/s)'); grid on;

nexttile; plot(t_sim, sim_history.T_applied/1e3, 'LineWidth', 1.3); 
ylabel('Thrust (kN)'); grid on;

nexttile; plot(t_sim, sim_history.m, 'LineWidth', 1.3); 
hold on; yline(P.m_dry, 'k--', 'LineWidth', 1);
ylabel('Mass (kg)'); grid on;

nexttile; plot(t_sim, sim_history.scp_iters, 'LineWidth', 1.3); 
ylabel('SCP iters'); xlabel('t (s)'); grid on;

linkaxes(fm.handles.MPC_Summary.Children.Children, 'x');

fm.exportAll;

%% ====================== Helper Functions ==========================================

function T_est = estimateTimeToTouchdown(state, P)
% Simple estimate of time to touchdown based on current altitude and velocity
% Assumes constant deceleration to zero velocity at zero altitude
p = state(1); v = state(2); % m = state(3);

if v >= 0  % moving up or stationary
    T_est = 20;  % arbitrary reasonable horizon
elseif p <= 0  % already on ground
    T_est = P.min_horizon_time;
else
    % Simple kinematic estimate: need to decelerate from v to 0 over distance p
    % Using average deceleration: a = v^2 / (2*p), time = v/a = 2*p/v
    if abs(v) > 1e-3
        T_est = 2 * p / abs(v);
        T_est = max(T_est, P.min_horizon_time);
        T_est = min(T_est, P.T_max_mission);  % cap at reasonable value
    else
        T_est = P.min_horizon_time;
    end
end
end

function [sol, log] = runSCP(initial_state, T_horizon, N, P, max_iters, dt_override)
% Run SCP optimization from given initial state over specified horizon

% Create temporary parameters for this SCP run
P_scp = P;
P_scp.p0 = initial_state(1);
P_scp.v0 = initial_state(2);
P_scp.m0 = initial_state(3);
P_scp.T = T_horizon;
P_scp.N = N;
if nargin >= 6 && ~isempty(dt_override)
    P_scp.dt = dt_override;
else
    P_scp.dt = T_horizon / N;
end
P_scp.max_iters = max_iters;  % use the adaptive max_iters

try
    % Initialize reference trajectory
    [p_ref, v_ref, T_ref, m_ref] = initialReference(P_scp);
    
    % Initialize logging
    log = struct('cost',[],'slack',[],'rho_T',[],'rho_v',[]);
    
    % Initialize trust regions
    rho_T = P_scp.trust_init_T * ones(P_scp.N,1);
    rho_v = P_scp.trust_init_v * ones(P_scp.N+1,1);
    
    % SCP iterations
    for iter = 1:P_scp.max_iters
        % Linearize dynamics
        lin = linearizeDynamics(v_ref, T_ref, m_ref, P_scp);
        
        % Build and solve subproblem
        prob = buildSubproblem(p_ref,v_ref,T_ref,m_ref,lin,rho_T,rho_v,P_scp,iter);
        [z, cost_val, exitflag] = solveQP(prob);
        
        if exitflag <= 0
            rho_T = rho_T * P_scp.trust_shrink; 
            rho_v = rho_v * P_scp.trust_shrink;
            continue;
        end
        
        sol = extractSolution(z, prob);
        
        % Logging
        slack_norm = sum(sol.s_plus + sol.s_minus);
        log.cost(end+1,1) = cost_val;
        log.slack(end+1,1) = slack_norm;
        log.rho_T(end+1,1) = mean(rho_T);
        log.rho_v(end+1,1) = mean(rho_v);
        
        % Check convergence
        if iter > 1
            rel_impr = abs(log.cost(end)-log.cost(end-1))/max(1,abs(log.cost(end-1)));
            if rel_impr < P_scp.tol_cost && slack_norm < P_scp.tol_slack
                break;
            end
        end
        
        % Trust region adaptation
        if slack_norm > P_scp.slack_trigger
            rho_T = rho_T * P_scp.trust_shrink; 
            rho_v = rho_v * P_scp.trust_shrink;
        else
            rho_T = rho_T * P_scp.trust_expand; 
            rho_v = rho_v * P_scp.trust_expand;
        end
        
        % Update reference
        p_ref = sol.p; v_ref = sol.v; T_ref = sol.T; m_ref = sol.m;
    end
    
catch ME
    fprintf('SCP error: %s\n', ME.message);
    sol = [];
    log = struct('cost',[],'slack',[],'rho_T',[],'rho_v',[]);
end

if isempty(sol)
    disp("")
end

end

function new_state = simulateOneStep(state, thrust, dt, P)
% Simulate one timestep forward using true nonlinear dynamics
p = state(1); v = state(2); m = state(3);

% Enforce thrust bounds and fuel constraints
T_k = min(max(thrust, P.T_min), P.T_max);
fuel_left = m - P.m_dry;

% Handle fuel depletion
if fuel_left <= 0
    T_k = 0;
    mdot = 0;
else
    mdot = T_k / (P.Isp * P.g0);
    if mdot * dt > fuel_left
        % Partial burn
        burn_frac = fuel_left / (mdot * dt);
        T_k = T_k * burn_frac;
        mdot = mdot * burn_frac;
    end
end

% Calculate forces
rho = isaDensity(p);
D = 0.5 * rho * P.Cd * P.A_ref * v * abs(v);  % drag force
g = localGravity(p, P);

% Equations of motion (using current mass)
a = (T_k - D) / m - g;

% Integrate using Euler method
v_new = v + dt * a;
p_new = p + dt * v;  % use current velocity for position update
m_new = max(P.m_dry, m - mdot * dt);

new_state = [p_new; v_new; m_new];
end

% Include all the helper functions from the original script
function [p_ref, v_ref, T_ref, m_ref] = initialReference(P)
% Construct a crude reference: thrust = mid between min & max until near ground, then taper.
T_ref = ones(P.N,1) * 0.6*(P.T_max + P.T_min)/2; % mid-high throttle
m_ref = zeros(P.N+1,1); m_ref(1)=P.m0;
for k=1:P.N
    mdot = T_ref(k)/(P.Isp*P.g0);
    m_ref(k+1) = max(P.m_dry, m_ref(k) - mdot*P.dt);
end
% Ballistic with average accel
p_ref = linspace(P.p0, P.p_target, P.N+1)';
v_ref = zeros(P.N+1,1); v_ref(1)=P.v0;
for k=1:P.N
    v_ref(k+1) = v_ref(k) + ( (P.v_target - P.v0)/P.T )*P.dt;
end
% Clip thrust to bounds explicitly (already mid-range)
T_ref = min(max(T_ref, P.T_min), P.T_max);
end

function lin = linearizeDynamics(v_ref, T_ref, m_ref, P)
% Linearize drag & mass flow about reference.
N = P.N; lin = struct();
% Drag: D = 0.5 * rho * Cd A v^2 sign(v) ; accel term a_d = -D / m
v_seg = v_ref(1:end-1);
rho_seg = isaDensity( (0.5*(v_ref(1:end-1)*0 + v_ref(2:end)*0) + 0) + 0 ); %#ok<NASGU>
% For density variation with altitude we need altitude; use p_ref? We linearize separately -> postpone to buildSubproblem using p_ref.
lin.v_ref = v_ref; lin.T_ref = T_ref; lin.m_ref = m_ref;
end

function prob = buildSubproblem(p_ref,v_ref,T_ref,m_ref,lin,rho_T,rho_v,P,iter)
% Variables: [p(0..N); v(0..N); m(0..N); T(0..N-1); s_plus(0..N-1); s_minus(0..N-1)]
N=P.N; n_p=N+1; n_v=N+1; n_m=N+1; n_T=N; n_sp=N; n_sm=N;
prob.idx.p  = 1:n_p;
prob.idx.v  = n_p + (1:n_v);
prob.idx.m  = n_p + n_v + (1:n_m);
prob.idx.T  = n_p + n_v + n_m + (1:n_T);
prob.idx.sp = n_p + n_v + n_m + n_T + (1:n_sp);
prob.idx.sm = n_p + n_v + n_m + n_T + n_sp + (1:n_sm);
prob.n_var  = prob.idx.sm(end);

% Preallocate cost
H = sparse(prob.n_var, prob.n_var); f = zeros(prob.n_var,1);
% Thrust effort (force^2)
H(prob.idx.T,prob.idx.T) = H(prob.idx.T,prob.idx.T) + 2*P.w_u*speye(n_T);
% Thrust rate smoothing
D = spdiags([-ones(n_T,1) [1; ones(n_T-1,1)]],[0 1],n_T-1,n_T); % forward diff
H(prob.idx.T,prob.idx.T) = H(prob.idx.T,prob.idx.T) + 2*P.w_du*(D'*D);
% Slack L1
f(prob.idx.sp) = P.w_s; f(prob.idx.sm) = P.w_s;

% Equality constraints count:
% Position dyn (N) + Velocity dyn (N) + Mass dyn (N) + initial p,v,m (3) + terminal p,v (2) = 3N +5
n_eq = 3*N + 5; Aeq = sparse(n_eq, prob.n_var); beq = zeros(n_eq,1); row=0;

for k=1:N
    % Position: p_{k+1} - p_k - dt v_k = 0
    row=row+1; Aeq(row,prob.idx.p(k))=-1; Aeq(row,prob.idx.p(k+1))=1; Aeq(row,prob.idx.v(k))=-P.dt;
end

for k=1:N
    % Velocity linearization: v_{k+1} = v_k + dt*( (T_k/m_k) - g(h_k) - D/m_k ) + slack
    % We'll approximate 1/m and drag terms affine around reference (m_ref, v_ref, p_ref).
    % Need altitude for density: use p_ref(k).
    h_k = p_ref(k);
    m_k_ref = m_ref(k);
    v_k_ref = v_ref(k);
    rho_k = isaDensity(h_k);
    D_ref = 0.5 * rho_k * P.Cd * P.A_ref * v_k_ref*abs(v_k_ref); % drag force magnitude
    a_thrust_coeff = P.dt / m_k_ref;      % derivative wrt T_k
    % Acceleration wrt v (drag linearization about v_ref): D ≈ D_ref + dD/dv * (v - v_ref)
    dD_dv = 0.5 * rho_k * P.Cd * P.A_ref * (2*abs(v_k_ref));
    drag_sign = sign(v_k_ref) + (v_k_ref==0); % to avoid zero sign -> treat as +1
    dD_dv = dD_dv * drag_sign; % directional derivative
    a_drag_const = D_ref / m_k_ref; % part of accel from drag
    a_drag_lin_coeff = dD_dv / m_k_ref; % multiplies (v_k - v_k_ref)
    g_k = localGravity(h_k, P);

    row=row+1;
    % LHS: v_{k+1} - v_k - dt*( T_k/m_ref - g_k - (a_drag_const + a_drag_lin_coeff*(v_k - v_k_ref)) ) - (s^+ - s^-) = 0
    Aeq(row,prob.idx.v(k))   = -1 - P.dt * (-a_drag_lin_coeff); % because -dt*(- a_drag_lin_coeff*v_k) = + dt * a_drag_lin_coeff
    Aeq(row,prob.idx.v(k+1)) =  1;
    Aeq(row,prob.idx.T(k))   = -a_thrust_coeff; % -dt*(T_k/m_ref) -> coefficient -dt/m_ref
    Aeq(row,prob.idx.sp(k))  = -1;
    Aeq(row,prob.idx.sm(k))  =  1;
    beq(row) = P.dt*( -g_k - a_drag_const + a_drag_lin_coeff*v_k_ref );
end

for k=1:N
    % Mass: m_{k+1} = m_k - dt * mdot_k ; mdot_k = T_k / (Isp*g0)
    row=row+1;
    Aeq(row,prob.idx.m(k))   = -1;
    Aeq(row,prob.idx.m(k+1)) =  1;
    Aeq(row,prob.idx.T(k))   =  P.dt/(P.Isp*P.g0); % since -dt * (T/(Isp g0)) -> move to LHS: + dt/(Isp g0) * T_k
end

% Initial conditions
row=row+1; Aeq(row,prob.idx.p(1))=1; beq(row)=P.p0;
row=row+1; Aeq(row,prob.idx.v(1))=1; beq(row)=P.v0;
row=row+1; Aeq(row,prob.idx.m(1))=1; beq(row)=P.m0;
% Terminal constraints
row=row+1; Aeq(row,prob.idx.p(end))=1; beq(row)=P.p_target;
row=row+1; Aeq(row,prob.idx.v(end))=1; beq(row)=P.v_target;

% Bounds
lb = -inf(prob.n_var,1); ub = inf(prob.n_var,1);
% Thrust bounds
lb(prob.idx.T) = P.T_min; ub(prob.idx.T) = P.T_max;
% Slack positivity
lb(prob.idx.sp) = 0; lb(prob.idx.sm) = 0;
% Mass lower bound
lb(prob.idx.m) = P.m_dry;

% Propellant exhaustion enforcement (approx): ensure total integrated mdot <= fuel.
% Soft way: after solve we check replay; hard way would add inequality sum(dt*T/(Isp g0)) <= fuel.
% Implement hard inequality via Aineq * z <= bineq.
Aineq = sparse(1, prob.n_var); bineq = P.fuel_mass;
for k=1:N
    Aineq(1, prob.idx.T(k)) = P.dt / (P.Isp * P.g0); % sum mdot*dt = fuel_used <= fuel_mass
end

% Trust regions (skip on first iter)
if iter>1
    ub(prob.idx.T) = min(ub(prob.idx.T), T_ref + rho_T);
    lb(prob.idx.T) = max(lb(prob.idx.T), T_ref - rho_T);
    ub(prob.idx.v) = min(ub(prob.idx.v), v_ref + rho_v);
    lb(prob.idx.v) = max(lb(prob.idx.v), v_ref - rho_v);
end

prob.H=H; prob.f=f; prob.Aeq=Aeq; prob.beq=beq; prob.lb=lb; prob.ub=ub; prob.Aineq=Aineq; prob.bineq=bineq; prob.P=P; prob.N=N;
end

function [z, cost_val, exitflag] = solveQP(prob)
opts = optimoptions('quadprog','Display','off','Algorithm','interior-point-convex');
[z, fval, exitflag] = quadprog(prob.H, prob.f, prob.Aineq, prob.bineq, prob.Aeq, prob.beq, prob.lb, prob.ub, [], opts);
cost_val = fval;
end

function sol = extractSolution(z, prob)
sol.p = z(prob.idx.p); sol.v = z(prob.idx.v); sol.m = z(prob.idx.m); sol.T = z(prob.idx.T); sol.s_plus = z(prob.idx.sp); sol.s_minus = z(prob.idx.sm);
end

function g = localGravity(h, P)
if P.use_var_g
    r = P.Re + h; g = P.mu / r^2; else; g = P.g0; end
end

function rho = isaDensity(h)
% Exponential model adequate below ~10 km for demo.
H_scale = 8500; rho0 = 1.225; rho = rho0 * exp(-h / H_scale);
end
