function create_scp_summary_plot(scp_debug_history, fm)
%CREATE_SCP_SUMMARY_PLOT Creates a summary plot of all SCP calls during simulation
%
% This function creates overview plots showing how SCP performance evolves
% throughout the entire MPC simulation, helping identify systematic issues.

if isempty(scp_debug_history.time)
    fprintf('Warning: No SCP debug history available\n');
    return;
end

n_scp_calls = length(scp_debug_history.time);
scp_times = scp_debug_history.time;

% Extract summary statistics from each SCP call
final_costs = zeros(n_scp_calls, 1);
final_slacks = zeros(n_scp_calls, 1);
n_iterations = zeros(n_scp_calls, 1);
mean_thrusts = zeros(n_scp_calls, 1);
mean_gimbals = zeros(n_scp_calls, 1);
cost_flatness = zeros(n_scp_calls, 1);  % Check if cost is flat
convergence_achieved = zeros(n_scp_calls, 1);

for i = 1:n_scp_calls
    log = scp_debug_history.log{i};
    sol = scp_debug_history.sol{i};
    
    final_costs(i) = log.cost(end);
    final_slacks(i) = log.slack(end);
    n_iterations(i) = length(log.cost);
    
    if isfield(log, 'thrust_mean')
        mean_thrusts(i) = log.thrust_mean(end);
        mean_gimbals(i) = log.gimbal_rms(end);
    end
    
    % Check if cost function is flat (major red flag)
    if length(log.cost) > 1
        cost_variation = std(log.cost) / max(abs(mean(log.cost)), 1e-12);
        cost_flatness(i) = cost_variation < 1e-10;  % Effectively flat
    end
    
    % Check if properly converged vs hit iteration limit
    if length(log.cost) > 1
        final_cost_change = abs(log.cost(end) - log.cost(end-1)) / max(abs(log.cost(end-1)), 1e-12);
        convergence_achieved(i) = (final_cost_change < 1e-6) && (log.slack(end) < 1e-6);
    end
end

%% Main Summary Figure
fh = fm.newFigure('SCP_Summary_All_Calls');
fig = fh.Parent;
TL = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(TL, sprintf('SCP Performance Summary (%d calls)', n_scp_calls));

% SCP Call Timeline
nexttile;
bar(1:n_scp_calls, scp_times, 'FaceColor', 'cyan');
xlabel('SCP Call Number'); ylabel('Simulation Time (s)');
title('SCP Call Timeline');
grid on;

% Final Cost Evolution
nexttile;
semilogy(1:n_scp_calls, final_costs, 'bo-', 'LineWidth', 2);
xlabel('SCP Call Number'); ylabel('Final Cost');
title('Final Cost per SCP Call');
grid on;

% Mark problematic calls
if any(cost_flatness)
    hold on;
    flat_calls = find(cost_flatness);
    semilogy(flat_calls, final_costs(flat_calls), 'rx', 'MarkerSize', 10, 'LineWidth', 3);
    legend('Normal', 'Flat Cost', 'Location', 'best');
end

% Slack Evolution
nexttile;
semilogy(1:n_scp_calls, final_slacks, 'ro-', 'LineWidth', 2);
xlabel('SCP Call Number'); ylabel('Final Slack');
title('Slack Variables per SCP Call');
grid on;

% Iterations per Call
nexttile;
bar(1:n_scp_calls, n_iterations, 'FaceColor', 'green');
xlabel('SCP Call Number'); ylabel('Iterations to Converge');
title('SCP Iteration Count');
grid on;

% Thrust Output Analysis
nexttile;
plot(1:n_scp_calls, mean_thrusts/1e3, 'bo-', 'LineWidth', 2);
xlabel('SCP Call Number'); ylabel('Mean Thrust (kN)');
title('Thrust Output per SCP Call');
grid on;

% Highlight zero thrust issues
if any(mean_thrusts < 1e3)
    hold on;
    zero_thrust_calls = find(mean_thrusts < 1e3);
    plot(zero_thrust_calls, mean_thrusts(zero_thrust_calls)/1e3, 'rx', 'MarkerSize', 10, 'LineWidth', 3);
    legend('Normal', 'Zero Thrust', 'Location', 'best');
end

% Gimbal Output Analysis
nexttile;
plot(1:n_scp_calls, rad2deg(mean_gimbals), 'go-', 'LineWidth', 2);
xlabel('SCP Call Number'); ylabel('RMS Gimbal (deg)');
title('Gimbal Output per SCP Call');
grid on;

% Problem Indicators
nexttile;
problem_score = cost_flatness + (mean_thrusts < 1e3) + (final_slacks > 1e-3);
bar(1:n_scp_calls, problem_score, 'FaceColor', 'red');
xlabel('SCP Call Number'); ylabel('Problem Score');
title('Problem Indicators (0=Good, 3=Bad)');
grid on;
ylim([0, 3.5]);

% Convergence Quality
nexttile;
bar(1:n_scp_calls, convergence_achieved, 'FaceColor', 'blue');
xlabel('SCP Call Number'); ylabel('Proper Convergence');
title('Convergence Quality (1=Good, 0=Poor)');
grid on;
ylim([0, 1.2]);

% Cost Function Health Check
nexttile;
bar(1:n_scp_calls, ~cost_flatness, 'FaceColor', 'magenta');
xlabel('SCP Call Number'); ylabel('Cost Function Active');
title('Cost Function Health (1=Active, 0=Flat)');
grid on;
ylim([0, 1.2]);

%% Diagnostic Summary Text

% Create diagnostic text
diagnostic_text = sprintf('DIAGNOSTIC SUMMARY:\\n\\n');

% Overall health check
healthy_calls = sum(~cost_flatness & mean_thrusts > 1e3 & final_slacks < 1e-3);
diagnostic_text = [diagnostic_text sprintf('Healthy SCP calls: %d/%d (%.1f%%)\\n', ...
    healthy_calls, n_scp_calls, 100*healthy_calls/n_scp_calls)];

% Specific issues
zero_thrust_count = sum(mean_thrusts < 1e3);
flat_cost_count = sum(cost_flatness);
high_slack_count = sum(final_slacks > 1e-3);

if zero_thrust_count > 0
    diagnostic_text = [diagnostic_text sprintf('*** ZERO THRUST ISSUE: %d calls ***\\n', zero_thrust_count)];
end
if flat_cost_count > 0
    diagnostic_text = [diagnostic_text sprintf('*** FLAT COST ISSUE: %d calls ***\\n', flat_cost_count)];
end
if high_slack_count > 0
    diagnostic_text = [diagnostic_text sprintf('*** HIGH SLACK ISSUE: %d calls ***\\n', high_slack_count)];
end

% Potential causes and recommendations
diagnostic_text = [diagnostic_text sprintf('\\nPOTENTIAL CAUSES:\\n')];
if flat_cost_count > 0.5 * n_scp_calls
    diagnostic_text = [diagnostic_text sprintf('- Cost function weights may be too small\\n')];
    diagnostic_text = [diagnostic_text sprintf('- Initial reference may be poor\\n')];
    diagnostic_text = [diagnostic_text sprintf('- Linearization may have errors\\n')];
end
if zero_thrust_count > 0.5 * n_scp_calls
    diagnostic_text = [diagnostic_text sprintf('- Terminal constraints may be unreachable\\n')];
    diagnostic_text = [diagnostic_text sprintf('- Trust regions may be too restrictive\\n')];
    diagnostic_text = [diagnostic_text sprintf('- Control penalty weights too high\\n')];
end

% Add as text annotation
text(0.05, 0.95, diagnostic_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 10, 'FontName', 'Courier', 'Interpreter', 'none');

fprintf('\n=== SCP SUMMARY DIAGNOSTICS ===\n');
fprintf('Total SCP calls: %d\n', n_scp_calls);
fprintf('Healthy calls: %d (%.1f%%)\n', healthy_calls, 100*healthy_calls/n_scp_calls);
fprintf('Zero thrust calls: %d\n', zero_thrust_count);
fprintf('Flat cost calls: %d\n', flat_cost_count);
fprintf('High slack calls: %d\n', high_slack_count);

end