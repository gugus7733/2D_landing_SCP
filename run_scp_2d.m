function [sol, log] = run_scp_2d(initial_state, T_horizon, N, P, dt_override, max_iters_override, prev_sol)
%RUN_SCP_2D Solves the 2D rocket landing problem using SCP.
%
% This function initializes a reference trajectory and then iteratively
% linearizes the dynamics and solves a convex subproblem (QP) to find an
% optimal state-control trajectory.

% SCP Parameters
P_scp = P;
P_scp.X0 = initial_state;
P_scp.T = T_horizon;
P_scp.N = N;

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
    % Pass both the previous solution and applied timestep information
    dt_applied = P.dt_scp;  % Timestep that was actually applied in simulation
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

% Initialize trust regions (adaptive from main script parameters)
trust_T = P.trust_init_T;
trust_delta = P.trust_init_delta;
trust_vx = P.trust_init_vx;
trust_vy = P.trust_init_vy;
trust_omega = P.trust_init_omega;

% SCP Iteration Loop
for iter = 1:P_scp.max_iters
    % Pass current trust regions to P_scp for this iteration
    P_scp.trust_T = trust_T;
    P_scp.trust_delta = trust_delta;
    P_scp.trust_vx = trust_vx;
    P_scp.trust_vy = trust_vy;
    P_scp.trust_omega = trust_omega;
    
    % Build and solve the convex subproblem (Quadratic Program)
    prob = build_subproblem_2d(X_ref, U_ref, P_scp);
    [z, cost_val, exitflag] = solve_qp(prob);

    if exitflag <= 0
        fprintf('  QP solver failed (iter %d)\n', iter);
        % Basic recovery: if solver fails, we keep the old reference and exit
        % A more advanced implementation would shrink trust regions here.
        if iter > 1
            sol = extract_solution(z_prev, prob);
        else
            sol = []; % Failed on first iteration
        end
        return;
    end

    % Extract solution from the solver's vector 'z'
    sol = extract_solution(z, prob);
    z_prev = z; % Store for recovery

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

    % Check for convergence
    if iter > 1
        cost_change = abs(log.cost(iter) - log.cost(iter-1)) / max(abs(log.cost(iter-1)), 1e-8);
        if cost_change < P.tol_cost && slack_norm < P.tol_slack
%             fprintf('    Converged: cost_change=%.2e, slack=%.2e\n', cost_change, slack_norm);
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

% Store final trust region values for inspection
log.final_trust_T = trust_T;
log.final_trust_delta = trust_delta;
log.final_trust_vx = trust_vx;
log.final_trust_vy = trust_vy;
log.final_trust_omega = trust_omega;

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