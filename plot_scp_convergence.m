function plot_scp_convergence(scp_log, scp_sol, sim_time, fm, P)
%PLOT_SCP_CONVERGENCE Creates detailed SCP convergence analysis plots
%
% This function creates comprehensive plots to debug SCP optimization issues:
% - Cost function evolution
% - Slack variable breakdown
% - Trust region adaptation
% - Solution quality metrics
% - Control output analysis

if isempty(scp_log) || isempty(scp_log.cost)
    fprintf('Warning: No SCP log data available for plotting\n');
    return;
end

n_iters = length(scp_log.cost);
iter_vec = 1:n_iters;

%% Main SCP Convergence Analysis Figure
fh = fm.newFigure('SCP_Convergence_Analysis');
fig = fh.Parent;
TL = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(TL, sprintf('SCP Convergence Analysis (t=%.2f s, %d iterations)', sim_time, n_iters));

% Cost Function Evolution
nexttile;
semilogy(iter_vec, scp_log.cost, 'b-o', 'LineWidth', 2);
grid on; ylabel('Cost Function'); xlabel('SCP Iteration');
title('Cost Evolution');
if n_iters > 1
    cost_change = abs(diff(scp_log.cost));
    if all(cost_change < 1e-12)
        text(0.1, 0.9, 'FLAT COST!', 'Units', 'normalized', 'Color', 'red', 'FontWeight', 'bold');
    end
end

% Slack Variables Breakdown
nexttile;
if isfield(scp_log, 'slack_vx')
    semilogy(iter_vec, scp_log.slack_vx, 'r-o', 'LineWidth', 1.5); hold on;
    semilogy(iter_vec, scp_log.slack_vy, 'g-o', 'LineWidth', 1.5);
    semilogy(iter_vec, scp_log.slack_omega, 'b-o', 'LineWidth', 1.5);
    legend('Slack vx', 'Slack vy', 'Slack \omega', 'Location', 'best');
else
    semilogy(iter_vec, scp_log.slack, 'k-o', 'LineWidth', 2);
end
grid on; ylabel('Slack Variables'); xlabel('SCP Iteration');
title('Slack Variable Evolution');

% Trust Region Evolution
nexttile;
if isfield(scp_log, 'trust_T_history')
    yyaxis left;
    plot(iter_vec, scp_log.trust_T_history/1e3, 'b-o', 'LineWidth', 2);
    ylabel('Trust T (kN)');
    yyaxis right;
    plot(iter_vec, rad2deg(scp_log.trust_delta_history), 'r-o', 'LineWidth', 2);
    ylabel('Trust \delta (deg)');
    xlabel('SCP Iteration');
end
grid on; title('Trust Region Evolution');

% Control Output Analysis
nexttile;
if isfield(scp_log, 'thrust_mean')
    plot(iter_vec, scp_log.thrust_mean/1e3, 'b-o', 'LineWidth', 2); hold on;
    plot(iter_vec, scp_log.thrust_max/1e3, 'r--o', 'LineWidth', 1.5);
    legend('Mean Thrust', 'Max Thrust', 'Location', 'best');
    ylabel('Thrust (kN)'); xlabel('SCP Iteration');
    title('Thrust Evolution');
    
    % Highlight zero thrust issue
    if all(scp_log.thrust_mean < 1e3)  % Less than 1 kN
        text(0.1, 0.9, 'ZERO THRUST!', 'Units', 'normalized', 'Color', 'red', 'FontWeight', 'bold');
    end
end
grid on;

% Gimbal Angle Analysis
nexttile;
if isfield(scp_log, 'gimbal_rms')
    plot(iter_vec, rad2deg(scp_log.gimbal_rms), 'g-o', 'LineWidth', 2);
    ylabel('Gimbal RMS (deg)'); xlabel('SCP Iteration');
    title('Gimbal Angle Evolution');
    
    % Highlight zero gimbal issue
    if all(rad2deg(scp_log.gimbal_rms) < 0.1)  % Less than 0.1 degrees
        text(0.1, 0.9, 'ZERO GIMBAL!', 'Units', 'normalized', 'Color', 'red', 'FontWeight', 'bold');
    end
end
grid on;

% Reference Change Analysis
nexttile;
if isfield(scp_log, 'X_ref_change')
    semilogy(iter_vec, max(scp_log.X_ref_change, 1e-16), 'b-o', 'LineWidth', 2); hold on;
    semilogy(iter_vec, max(scp_log.U_ref_change, 1e-16), 'r-o', 'LineWidth', 2);
    legend('State Change', 'Control Change', 'Location', 'best');
    ylabel('Reference Update Norm'); xlabel('SCP Iteration');
    title('Reference Update Magnitude');
    
    % Check for stagnation
    if n_iters > 2 && all(scp_log.X_ref_change(3:end) < 1e-12)
        text(0.1, 0.9, 'NO UPDATES!', 'Units', 'normalized', 'Color', 'red', 'FontWeight', 'bold');
    end
end
grid on;

% QP Solver Status
nexttile;
if isfield(scp_log, 'qp_exitflag')
    bar(iter_vec, scp_log.qp_exitflag, 'FaceColor', 'cyan');
    ylabel('QP Exit Flag'); xlabel('SCP Iteration');
    title('QP Solver Status');
    ylim([-1, 4]);
    grid on;
    
    % Add status interpretation
    if any(scp_log.qp_exitflag <= 0)
        text(0.1, 0.9, 'QP FAILURES!', 'Units', 'normalized', 'Color', 'red', 'FontWeight', 'bold');
    end
end

% Cost Change Analysis
nexttile;
if n_iters > 1
    cost_change = [0, abs(diff(scp_log.cost))];
    semilogy(iter_vec, max(cost_change, 1e-16), 'k-o', 'LineWidth', 2);
    ylabel('|Cost Change|'); xlabel('SCP Iteration');
    title('Cost Change Magnitude');
    grid on;
    
    % Highlight convergence threshold
    if isfield(scp_log, 'tol_cost')
        yline(scp_log.tol_cost, 'r--', 'Convergence Threshold');
    end
end

% Convergence Rate Analysis
nexttile;
if n_iters > 2
    % Compute relative improvement
    rel_impr = zeros(size(iter_vec));
    for i = 2:n_iters
        if abs(scp_log.cost(i-1)) > 1e-12
            rel_impr(i) = abs(scp_log.cost(i) - scp_log.cost(i-1)) / abs(scp_log.cost(i-1));
        end
    end
    semilogy(iter_vec, max(rel_impr, 1e-16), 'm-o', 'LineWidth', 2);
    ylabel('Relative Cost Change'); xlabel('SCP Iteration');
    title('Convergence Rate');
    grid on;
end

%% Solution Quality Analysis (if solution exists)
if ~isempty(scp_sol)
    fh2 = fm.newFigure('SCP_Solution_Analysis');
    fig2 = fh2.Parent;
    TL2 = tiledlayout(fig2, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(TL2, sprintf('SCP Solution Quality (t=%.2f s)', sim_time));
    
    N = size(scp_sol.U, 2);
    % Use actual SCP timestep and compute horizon time
    dt_scp = P.dt_scp;
    T_horizon = N * dt_scp;
    time_vec = linspace(0, T_horizon, N);
    
    fprintf('  SCP Solution: N=%d steps, dt=%.2fs, horizon=%.2fs\n', N, dt_scp, T_horizon);
    
    % Control Trajectories
    nexttile;
    plot(time_vec, scp_sol.U(1,:)/1e3, 'b-', 'LineWidth', 2); hold on;
    ylabel('Thrust (kN)'); xlabel('Time in Horizon (s)');
    title('Optimal Thrust Profile');
    grid on;
    
    nexttile;
    plot(time_vec, rad2deg(scp_sol.U(2,:)), 'r-', 'LineWidth', 2);
    ylabel('Gimbal Angle (deg)'); xlabel('Time in Horizon (s)');
    title('Optimal Gimbal Profile');
    grid on;
    
    % State Trajectories (key states)
    nexttile;
    plot(time_vec, scp_sol.X(1,1:end-1), 'b-', 'LineWidth', 2); hold on;
    plot(time_vec, scp_sol.X(2,1:end-1), 'r-', 'LineWidth', 2);
    legend('x (m)', 'y (m)', 'Location', 'best');
    ylabel('Position'); xlabel('Time in Horizon (s)');
    title('Predicted Position');
    grid on;
    
    nexttile;
    plot(time_vec, scp_sol.X(3,1:end-1), 'b-', 'LineWidth', 2); hold on;
    plot(time_vec, scp_sol.X(4,1:end-1), 'r-', 'LineWidth', 2);
    legend('vx (m/s)', 'vy (m/s)', 'Location', 'best');
    ylabel('Velocity'); xlabel('Time in Horizon (s)');
    title('Predicted Velocity');
    grid on;
end

% Diagnostic summary removed per user request

end