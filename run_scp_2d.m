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

log = struct('cost',[], 'slack',[], 'retry_level', [], 'thrust_mean',[], 'thrust_max',[], 'gimbal_rms',[], ...
            'X_ref_change',[], 'U_ref_change',[], 'trust_T_history',[], 'trust_delta_history',[], ...
            'trust_vx_history',[], 'trust_vy_history',[], 'trust_omega_history',[], ...
            'slack_vx',[], 'slack_vy',[], 'slack_omega',[], 'qp_exitflag',[], ...
            'merit_pred',[], 'merit_act',[], 'merit_ratio',[], 'trust_update_factor',[]);
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

% Initialize merit tracking variables
prev_merit = []; % Previous iteration's actual merit (empty for first iteration)

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
    [prob, z, cost_val, exitflag, retry_level] = solve_qp_with_retry(X_ref, U_ref, P_scp, t_remaining, iter);

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
    sol = extract_solution(z, prob);
    z_prev = z; % Store for recovery
    prob_prev = prob; % Store problem for recovery

    pred_merit = cost_val;
    
    % Compute actual merit at the new solution
    act_merit = compute_merit(sol.X, sol.U, sol.s_v, sol.s_w, ...
                             sol.s_T_upper, sol.s_T_lower, sol.s_delta_upper, sol.s_delta_lower, sol.s_terminal, ...
                             P_scp, t_remaining, iter);
    
    % Compute merit ratio for trust region update
    if ~isempty(prev_merit) % Skip first iteration

        rho_num = prev_merit - act_merit;      % actual reduction
        rho_den = pred_merit - prev_merit; % predicted reduction
        rho     = rho_num / max(abs(rho_den),1e-12);
        merit_ratio = rho;

        % Trust region update based on merit ratio ρ
        if merit_ratio < 0.25
            trust_update_factor = 0.5; % Shrink
            update_desc = 'shrink';
        elseif merit_ratio > 0.90
            trust_update_factor = 1.5; % Expand
            update_desc = 'expand';
        else
            trust_update_factor = 1.0; % Maintain
            update_desc = 'maintain';
        end
        
        % Apply trust region updates with clamping and safety checks
        
        trust_T = max(P.trust_min_T, min(P.trust_max_T, trust_T * trust_update_factor));
        trust_delta = max(P.trust_min_delta, min(P.trust_max_delta, trust_delta * trust_update_factor));
        trust_vx = max(P.trust_min_vx, min(P.trust_max_vx, trust_vx * trust_update_factor));
        trust_vy = max(P.trust_min_vy, min(P.trust_max_vy, trust_vy * trust_update_factor));
        trust_omega = max(P.trust_min_omega, min(P.trust_max_omega, trust_omega * trust_update_factor));
        
        % Log trust region changes for debugging
%         fprintf('    Trust region changes: T: %.2e→%.2e, δ: %.2e→%.2e, vₓ: %.1f→%.1f, vᵧ: %.1f→%.1f, ω: %.3f→%.3f\n', ...
%                 old_trust_T, trust_T, old_trust_delta, trust_delta, ...
%                 old_trust_vx, trust_vx, old_trust_vy, trust_vy, ...
%                 old_trust_omega, trust_omega);
%         
%             fprintf('  Merit ratio ρ=%.3f → %s trust regions (factor=%.1f)\n', merit_ratio, update_desc, trust_update_factor);
    else
        merit_ratio = NaN; % Invalid for first iteration
        trust_update_factor = 1.0;
%         fprintf('  First iteration: keeping initial trust regions\n');
    end

    % Log progress - zero out constraint violation slacks for fair merit calculation
    slack_norm = sum(abs(sol.s_v(:))) + sum(abs(sol.s_w(:)));
    % Note: s_T_upper, s_T_lower, s_delta_upper, s_delta_lower, s_terminal are not included
    % in merit calculation as they are temporary aids for early iterations
    log.cost(iter) = cost_val;
    log.slack(iter) = slack_norm;
    log.retry_level(iter) = retry_level;
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
    
    % Log merit function tracking
    log.merit_pred(iter) = pred_merit;
    log.merit_act(iter) = act_merit;
    log.merit_ratio(iter) = merit_ratio;
    log.trust_update_factor(iter) = trust_update_factor;

    log.final_trust_T = 0;
    log.final_trust_delta = 0;
    log.final_trust_vx = 0;
    log.final_trust_vy = 0;
    log.final_trust_omega = 0;
    
    % Update previous merit for next iteration
    prev_merit = act_merit;
    
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

    % Store trust regions in solution for persistence between SCP calls
    sol.trust_regions = struct('T', trust_T, 'delta', trust_delta, ...
                              'vx', trust_vx, 'vy', trust_vy, 'omega', trust_omega);

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

end


function [prob, z, cost_val, exitflag, retry_level] = solve_qp_with_retry(X_ref, U_ref, P_scp, t_remaining, iter)
    % Enhanced QP solver with automatic retry, slack relaxation, and trust region widening
    
    retry_level = 0;
    exitflag = -1; % Initialize as failed
    
    % Store original trust radii for restoration after retry
    orig_trust_T = P_scp.trust_T;
    orig_trust_delta = P_scp.trust_delta;
    orig_trust_vx = P_scp.trust_vx;
    orig_trust_vy = P_scp.trust_vy;
    orig_trust_omega = P_scp.trust_omega;
    
    for retry = 0:P_scp.max_slack_retries
        retry_level = retry;
        
        % Progressive trust region widening for retries
        if retry > 0
            retry_scale = 2^retry;
            
            % Widen trust regions with clamping (use safeguards for missing max parameters)
            if isfield(P_scp, 'trust_max_T')
                max_T = P_scp.trust_max_T;
            else
                max_T = 2*P_scp.T_max;
            end
            if isfield(P_scp, 'trust_max_delta')
                max_delta = P_scp.trust_max_delta;
            else
                max_delta = P_scp.delta_max;
            end
            if isfield(P_scp, 'trust_max_vx')
                max_vx = P_scp.trust_max_vx;
            else
                max_vx = 500;
            end
            if isfield(P_scp, 'trust_max_vy')
                max_vy = P_scp.trust_max_vy;
            else
                max_vy = 500;
            end
            if isfield(P_scp, 'trust_max_omega')
                max_omega = P_scp.trust_max_omega;
            else
                max_omega = deg2rad(180);
            end
            
            P_scp.trust_T = min(orig_trust_T * retry_scale, max_T);
            P_scp.trust_delta = min(orig_trust_delta * retry_scale, max_delta);
            P_scp.trust_vx = min(orig_trust_vx * retry_scale, max_vx);
            P_scp.trust_vy = min(orig_trust_vy * retry_scale, max_vy);
            P_scp.trust_omega = min(orig_trust_omega * retry_scale, max_omega);
            
%             fprintf('    Retry %d: widening trust regions by factor %.1f\n', retry, retry_scale);
%             fprintf('      T: %.1e → %.1e, δ: %.1e → %.1e, vₓ: %.1e → %.1e, vᵧ: %.1e → %.1e, ω: %.1e → %.1e\n', ...
%                     orig_trust_T, P_scp.trust_T, orig_trust_delta, P_scp.trust_delta, ...
%                     orig_trust_vx, P_scp.trust_vx, orig_trust_vy, P_scp.trust_vy, ...
%                     orig_trust_omega, P_scp.trust_omega);
        end
        
        % Build subproblem with current retry level (affects slack bounds) and possibly widened trust radii
        prob = build_subproblem_2d(X_ref, U_ref, P_scp, t_remaining, retry, iter);
        
        % Solve QP
        [z_attempt, cost_val, exitflag] = solve_qp(prob);
        
        if exitflag > 0
            % Success - return result
            z = z_attempt;
            if retry > 0
%                 fprintf('    QP succeeded on retry %d with relaxed constraints\n', retry);
            end
            
            % Restore original trust radii (merit-based updates will handle changes)
            P_scp.trust_T = orig_trust_T;
            P_scp.trust_delta = orig_trust_delta;
            P_scp.trust_vx = orig_trust_vx;
            P_scp.trust_vy = orig_trust_vy;
            P_scp.trust_omega = orig_trust_omega;
            
            return;
        else
            % Failure - log and prepare for retry
            if retry < P_scp.max_slack_retries
%                 fprintf('    QP failed (exitflag=%d), retrying with further relaxation (%d/%d)\n', ...
%                         exitflag, retry+1, P_scp.max_slack_retries);
            else
%                 fprintf('    QP failed (exitflag=%d) - all retries exhausted\n', exitflag);
            end
        end
    end
    
    % If we reach here, all retries failed - restore original trust radii
    P_scp.trust_T = orig_trust_T;
    P_scp.trust_delta = orig_trust_delta;
    P_scp.trust_vx = orig_trust_vx;
    P_scp.trust_vy = orig_trust_vy;
    P_scp.trust_omega = orig_trust_omega;
    
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
    
    % Extract constraint violation slacks (for softening)
    sol.s_T_upper = z(prob.idx.s_T_upper)';
    sol.s_T_lower = z(prob.idx.s_T_lower)';
    sol.s_delta_upper = z(prob.idx.s_delta_upper)';
    sol.s_delta_lower = z(prob.idx.s_delta_lower)';
    sol.s_terminal = z(prob.idx.s_terminal);
    
    sol.P_scp = P;

    sol.trust_regions = struct();
end

function merit = compute_merit(X, U, s_v, s_w, s_T_upper, s_T_lower, s_delta_upper, s_delta_lower, s_terminal, P, t_remaining, iter)
%COMPUTE_MERIT Calculate merit function L = J_quad + w·‖slack‖₁
%
% The merit function includes:
% 1. Quadratic cost from controls (fuel + effort + rate penalties)
% 2. Dynamics slack penalties (time-dependent weights)
% 3. Constraint violation slack penalties (large fixed weights)
%
% INPUTS:
%   X, U        - State and control trajectories
%   s_*         - Various slack variables
%   P           - Parameter structure
%   t_remaining - Time remaining for slack weight calculation
%   iter        - Current iteration (affects constraint violation penalty)

    N = P.N;
    dt = P.dt;
    
    % --- Quadratic Cost Component ---
    J_quad = 0;
    
    % Control magnitude penalties
    J_quad = J_quad + P.w_T * sum(U(1,:).^2);
    J_quad = J_quad + P.w_delta * sum(U(2,:).^2);
    
    % Control rate penalties
    if N > 1
        dT_dt = diff(U(1,:)) / dt;
        ddelta_dt = diff(U(2,:)) / dt;
        J_quad = J_quad + P.w_dT * sum(dT_dt.^2);
        J_quad = J_quad + P.w_ddelta * sum(ddelta_dt.^2);
    end
    
    % Angular velocity magnitude penalty
    if isfield(P, 'w_omega') && P.w_omega > 0
        omega_trajectory = X(6, :);  % ω is 6th state
        J_quad = J_quad + P.w_omega * sum(omega_trajectory.^2);
    end
    
    % --- Dynamics Slack Penalty Component ---
    % Get time-dependent slack weight
    w_slack_current = slack_management_utils('compute_slack_weight', t_remaining, P.T, P);
    
    dynamics_slack_penalty = w_slack_current * (sum(abs(s_v(:))) + sum(abs(s_w(:))));
    
    % --- Constraint Violation Slack Penalty Component ---
    constraint_slack_penalty = 1e6; % Large penalty to discourage use
    
    if iter <= 2 % Only active during softening iterations
        constraint_violation_penalty = constraint_slack_penalty * (...
            sum(abs(s_T_upper)) + sum(abs(s_T_lower)) + ...
            sum(abs(s_delta_upper)) + sum(abs(s_delta_lower)) + ...
            sum(abs(s_terminal)));
    else
        constraint_violation_penalty = 0; % No constraint violation slacks in later iterations
    end
    
    % --- Total Merit Function ---
    merit = J_quad + dynamics_slack_penalty + constraint_violation_penalty;
end