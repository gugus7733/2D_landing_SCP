function [sol, log] = run_scp_2d(initial_state, T_horizon, N, P, dt_override, max_iters_override, prev_sol, sim_time)
%RUN_SCP_2D Solves the 2D rocket landing problem using SCP.
%
% This function initializes a reference trajectory and then iteratively
% linearizes the dynamics and solves a convex subproblem (QP) to find an
% optimal state-control trajectory.
%
% INPUTS:
%   sim_time    - Current simulation time for accurate slack management (optional)

% SCP Parameters
P_scp = P;
P_scp.X0 = initial_state;
P_scp.T = T_horizon;
P_scp.N = N;

% Validate slack management parameters
[is_valid, error_msg] = slack_management_utils('validate_parameters', P);
if ~is_valid
    error('Slack management parameter validation failed: %s', error_msg);
end

% Use override parameters if provided, otherwise use defaults
if nargin >= 5 && ~isempty(dt_override)
    P_scp.dt = dt_override;
else
    P_scp.dt = T_horizon / N;
end

if nargin >= 6 && ~isempty(max_iters_override)
    P_scp.max_iters = max_iters_override;
else
    P_scp.max_iters = P.max_iters_scp;
end

% State and Control dimensions
P_scp.n_states = 7;
P_scp.n_controls = 2;

% Initialize reference trajectory with warm start
if nargin >= 7 && ~isempty(prev_sol)
    % Warm start: use previous solution shifted forward in time
    % Extract actual applied timestep from previous solution (more accurate)
    if isfield(prev_sol, 'P_scp') && isfield(prev_sol.P_scp, 'dt')
        dt_applied = prev_sol.P_scp.dt;
    else
        dt_applied = P.dt_scp; % Fallback
    end
    [X_ref, U_ref] = get_initial_reference_2d(P_scp, prev_sol, dt_applied);
else
    % Cold start: generate initial reference from scratch
    [X_ref, U_ref] = get_initial_reference_2d(P_scp);
end

log = struct('cost',[], 'slack',[], 'thrust_mean',[], 'thrust_max',[], 'gimbal_rms',[], ...
            'X_ref_change',[], 'U_ref_change',[], 'trust_T_history',[], 'trust_delta_history',[], ...
            'trust_vx_history',[], 'trust_vy_history',[], 'trust_omega_history',[], ...
            'slack_vx',[], 'slack_vy',[], 'slack_omega',[], 'qp_exitflag',[]);
sol = [];

% Initialize trust regions with persistence from previous solution
if nargin >= 7 && ~isempty(prev_sol) && isfield(prev_sol, 'trust_regions')
    % Use trust regions from previous solution (warm start)
    trust_T = prev_sol.trust_regions.T;
    trust_delta = prev_sol.trust_regions.delta;
    trust_vx = prev_sol.trust_regions.vx;
    trust_vy = prev_sol.trust_regions.vy;
    trust_omega = prev_sol.trust_regions.omega;
else
    % Initialize trust regions from parameters (cold start)
    trust_T = P.trust_init_T;
    trust_delta = P.trust_init_delta;
    trust_vx = P.trust_init_vx;
    trust_vy = P.trust_init_vy;
    trust_omega = P.trust_init_omega;
end

% SCP Iteration Loop with Enhanced Slack Management
for iter = 1:P_scp.max_iters
    % Pass current trust regions to P_scp for this iteration
    P_scp.trust_T = trust_T;
    P_scp.trust_delta = trust_delta;
    P_scp.trust_vx = trust_vx;
    P_scp.trust_vy = trust_vy;
    P_scp.trust_omega = trust_omega;
    
    % Calculate remaining time for slack management
    if nargin >= 8 && ~isempty(sim_time)
        % Use actual simulation time for accurate calculation
        if isfield(P, 'T_max_mission')
            t_remaining = P.T_max_mission - sim_time;
        else
            % Estimate time to touchdown from current state
            altitude = initial_state(2);
            vy = initial_state(4);
            if vy < -1.0 && altitude > 0
                t_est_touchdown = altitude / abs(vy);
                t_remaining = min(T_horizon, t_est_touchdown);
            else
                t_remaining = T_horizon;
            end
        end
        t_remaining = max(t_remaining, 0.1); % Ensure positive value
    else
        % Fallback to horizon-based estimation
        if iter == 1
            t_remaining = T_horizon;
        else
            t_remaining = T_horizon * 0.8; % More conservative than before
        end
    end
    
    % Enhanced QP solving with retry mechanism
    [z, cost_val, exitflag, retry_level] = solve_qp_with_retry(X_ref, U_ref, P_scp, t_remaining);

    if exitflag <= 0
        fprintf('  QP solver failed after %d retries (iter %d)\n', P.max_slack_retries, iter);
        % Advanced recovery: if all retries fail, we keep the old reference and exit
        if iter > 1
            sol = extract_solution(z_prev, prob_prev);
        else
            sol = []; % Failed on first iteration - no previous solution
        end
        return;
    end

    % Extract solution from the solver's vector 'z' 
    % Need to rebuild prob for extraction since it may have been modified in retries
    prob = build_subproblem_2d(X_ref, U_ref, P_scp, t_remaining, 0);
    sol = extract_solution(z, prob);
    z_prev = z; % Store for recovery
    prob_prev = prob; % Store problem for recovery

    % Log progress
    slack_norm = sum(abs(sol.s_v(:))) + sum(abs(sol.s_w(:)));
    log.cost(iter) = cost_val;
    log.slack(iter) = slack_norm;
    log.thrust_mean(iter) = mean(sol.U(1,:));
    log.thrust_max(iter) = max(sol.U(1,:));
    log.gimbal_rms(iter) = sqrt(mean(sol.U(2,:).^2));
    log.slack_vx(iter) = sum(abs(sol.s_v(1,:)));
    log.slack_vy(iter) = sum(abs(sol.s_v(2,:)));
    log.slack_omega(iter) = sum(abs(sol.s_w(:)));
    log.trust_T_history(iter) = trust_T;
    log.trust_delta_history(iter) = trust_delta;
    log.trust_vx_history(iter) = trust_vx;
    log.trust_vy_history(iter) = trust_vy;
    log.trust_omega_history(iter) = trust_omega;
    log.qp_exitflag(iter) = exitflag;
    
    % Debug cost balance (first iteration only)
    if iter == 1
%         fprintf('    Cost components: Thrust=%.1f kN, Gimbal=%.2f deg, Slack=%.2e\n', ...
%                 log.thrust_mean(iter)/1e3, rad2deg(log.gimbal_rms(iter)), slack_norm);
    end
    
    if iter > 1
        log.X_ref_change(iter) = norm(sol.X(:) - X_ref(:));
        log.U_ref_change(iter) = norm(sol.U(:) - U_ref(:));
    else
        log.X_ref_change(iter) = 0;
        log.U_ref_change(iter) = 0;
    end

    % --- Trust region adaptation (matching 1D model logic) ---
    if slack_norm > P.slack_trigger
        % High slack: shrink trust regions (reduce freedom, match 1D logic)
        trust_T = max(P.trust_min_T, trust_T * P.trust_shrink);
        trust_delta = max(P.trust_min_delta, trust_delta * P.trust_shrink);
        trust_vx = max(P.trust_min_vx, trust_vx * P.trust_shrink);
        trust_vy = max(P.trust_min_vy, trust_vy * P.trust_shrink);
        trust_omega = max(P.trust_min_omega, trust_omega * P.trust_shrink);
%         fprintf('    Slack high (%.2e), shrinking trust: T=%.2f kN, delta=%.2f deg, vx=%.1f m/s, vy=%.1f m/s, w=%.2f deg/s\n', ...
%                 slack_norm, trust_T/1e3, rad2deg(trust_delta), trust_vx, trust_vy, rad2deg(trust_omega));
    else
        % Low slack: expand trust regions (allow more freedom)
        trust_T = min(P.trust_max_T, trust_T * P.trust_expand);
        trust_delta = min(P.trust_max_delta, trust_delta * P.trust_expand);
        trust_vx = min(P.trust_max_vx, trust_vx * P.trust_expand);
        trust_vy = min(P.trust_max_vy, trust_vy * P.trust_expand);
        trust_omega = min(P.trust_max_omega, trust_omega * P.trust_expand);
%         fprintf('    Slack low (%.2e), expanding trust: T=%.2f kN, delta=%.2f deg, vx=%.1f m/s, vy=%.1f m/s, w=%.2f deg/s\n', ...
%                 slack_norm, trust_T/1e3, rad2deg(trust_delta), trust_vx, trust_vy, rad2deg(trust_omega));
    end

    % Check for convergence with time-dependent slack tolerance
    if iter > 1
        cost_change = abs(log.cost(iter) - log.cost(iter-1)) / max(abs(log.cost(iter-1)), 1e-8);
        tol_slack_current = slack_management_utils('compute_slack_tolerance', t_remaining, T_horizon, P);
        
        if cost_change < P.tol_cost && slack_norm < tol_slack_current
%             fprintf('    Converged: cost_change=%.2e, slack=%.2e (tol=%.2e)\n', ...
%                     cost_change, slack_norm, tol_slack_current);
            break; % Converged
        end
    end

    % Update reference trajectory for the next iteration
    X_ref = sol.X;
    U_ref = sol.U;
end

if iter == P_scp.max_iters
%     fprintf('  Warning: SCP reached max iterations (%d).\n', P_scp.max_iters);
end

% Store final trust region values in log for inspection
log.final_trust_T = trust_T;
log.final_trust_delta = trust_delta;
log.final_trust_vx = trust_vx;
log.final_trust_vy = trust_vy;
log.final_trust_omega = trust_omega;

% Store trust regions in solution for persistence between SCP calls
if ~isempty(sol)
    sol.trust_regions = struct('T', trust_T, 'delta', trust_delta, ...
                              'vx', trust_vx, 'vy', trust_vy, 'omega', trust_omega);
end

end


function [z, cost_val, exitflag, retry_level] = solve_qp_with_retry(X_ref, U_ref, P_scp, t_remaining)
    % Enhanced QP solver with automatic retry and slack relaxation
    
    retry_level = 0;
    exitflag = -1; % Initialize as failed
    
    for retry = 0:P_scp.max_slack_retries
        retry_level = retry;
        
        % Build subproblem with current retry level (affects slack bounds)
        prob = build_subproblem_2d(X_ref, U_ref, P_scp, t_remaining, retry);
        
        % Solve QP
        [z_attempt, cost_val, exitflag] = solve_qp(prob);
        
        if exitflag > 0
            % Success - return result
            z = z_attempt;
            if retry > 0
                fprintf('    QP succeeded on retry %d with relaxed slack bounds\n', retry);
            end
            return;
        else
            % Failure - log and prepare for retry
            if retry < P_scp.max_slack_retries
                fprintf('    QP failed (exitflag=%d), retrying with relaxed slack bounds (%d/%d)\n', ...
                        exitflag, retry+1, P_scp.max_slack_retries);
            else
                fprintf('    QP failed (exitflag=%d) - all retries exhausted\n', exitflag);
            end
        end
    end
    
    % If we reach here, all retries failed
    z = [];
    cost_val = inf;
end

function [z, cost_val, exitflag] = solve_qp(prob)
    opts = optimoptions('quadprog','Display','off','Algorithm','interior-point-convex');
    [z, fval, exitflag] = quadprog(prob.H, prob.f, prob.Aineq, prob.bineq, prob.Aeq, prob.beq, prob.lb, prob.ub, [], opts);
    cost_val = fval;
end

function sol = extract_solution(z, prob)
    P = prob.P;
    sol.X = reshape(z(prob.idx.X), P.n_states, P.N + 1);
    sol.U = reshape(z(prob.idx.U), P.n_controls, P.N);
    sol.s_v = reshape(z(prob.idx.s_v), 2, P.N); % [vx; vy] slacks
    sol.s_w = reshape(z(prob.idx.s_w), 1, P.N);
    sol.P_scp = P;
end