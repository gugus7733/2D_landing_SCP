function plot_open_loop_comparison(predicted_traj, achieved_traj, start_time, dt_sequence, P, fm)
%PLOT_OPEN_LOOP_COMPARISON Plots predicted vs achieved trajectories in open-loop mode
%
% INPUTS:
%   predicted_traj  - SCP predicted state trajectory [7 x N+1]
%   achieved_traj   - Achieved trajectory structure with .t and .X fields
%   start_time      - Open-loop start time (s)
%   dt_sequence     - Timestep sequence used for predictions
%   P               - Parameters structure
%   fm              - Figure manager

if isempty(predicted_traj) || isempty(achieved_traj.X)
    fprintf('Warning: No trajectory data available for open-loop comparison plot\n');
    return;
end

% Create time vectors for predicted trajectory
n_pred_points = size(predicted_traj, 2);
t_pred = start_time + (0:n_pred_points-1) * dt_sequence(1);

% Extract achieved trajectory data
t_achieved = achieved_traj.t;
X_achieved = achieved_traj.X;

% State labels and units
state_labels = {'X Position', 'Y Position', 'X Velocity', 'Y Velocity', 'Pitch Angle', 'Angular Velocity'};
state_units = {'(m)', '(m)', '(m/s)', '(m/s)', '(deg)', '(deg/s)'};
state_indices = [1, 2, 3, 4, 5, 6]; % Skip mass (index 7)

% Create comparison plot
fig = fm.newFigure('OpenLoop_Trajectory_Comparison');

% Set up 3x2 subplot layout
for i = 1:6
    subplot(3, 2, i);
    
    % Extract state data
    state_idx = state_indices(i);
    pred_data = predicted_traj(state_idx, :);
    achieved_data = X_achieved(state_idx, :);
    
    % Convert angles to degrees for display
    if state_idx == 5 || state_idx == 6
        pred_data = rad2deg(pred_data);
        achieved_data = rad2deg(achieved_data);
    end
    
    % Plot predicted trajectory
    plot(t_pred, pred_data, 'b-', 'LineWidth', 2.5, 'DisplayName', 'SCP Predicted');
    hold on;
    
    % Plot achieved trajectory
    plot(t_achieved, achieved_data, 'r--', 'LineWidth', 2, 'DisplayName', 'Achieved');
    
    % Mark open-loop start time
    yl = ylim;
    plot([start_time, start_time], yl, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Open-Loop Start');
    
    % Formatting
    grid on;
    xlabel('Time (s)');
    ylabel([state_labels{i} ' ' state_units{i}]);
    title([state_labels{i} ' Comparison'], 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    
    % Set consistent time limits for all subplots
    xlim([min([t_pred(1), t_achieved(1)]), max([t_pred(end), t_achieved(end)])]);
    
    % Add trajectory endpoint markers
    plot(t_pred(end), pred_data(end), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
    plot(t_achieved(end), achieved_data(end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
end

% Add overall title
sgtitle('Open-Loop Trajectory Verification: SCP Prediction vs Reality', ...
        'FontSize', 16, 'FontWeight', 'bold');

end