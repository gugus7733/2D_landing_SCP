function plot_interactive_trajectory_2d(sim_history, P, fm)
%PLOT_INTERACTIVE_TRAJECTORY_2D Interactive 2D rocket trajectory visualization
%
% INPUTS:
%   sim_history - Structure containing trajectory data with fields:
%                 .t - time vector
%                 .X - state matrix [x, y, vx, vy, theta, omega, mass]
%                 .U - control matrix [T, delta]
%                 .vectors - structure with .velocity, .orientation, .thrust
%   P           - Parameter structure
%   fm          - FigureManager instance
%
% OUTPUTS:
%   Creates interactive plot with playback controls

%% Extract trajectory data
t_data        = sim_history.t;
X_data        = sim_history.X;
velocity_data = sim_history.vectors.velocity;
orient_data   = sim_history.vectors.orientation;
thrust_data   = sim_history.vectors.thrust;

n_points = length(t_data);
if n_points < 2
    fprintf('Warning: Insufficient trajectory data for interactive plot\n');
    return;
end

%% Create figure using FigureManager
fig_handle = fm.newFigure('Interactive_Trajectory_2D');
% set(fig_handle, 'Position', [100, 100, 1200, 800]);
plotedit(fig_handle,'off')
clf(fig_handle);

%% Setup main axes
ax_main = axes('Parent', fig_handle, 'Position', [0.1, 0.2, 0.8, 0.7]);
hold(ax_main, 'on');
grid(ax_main, 'on');
axis(ax_main, 'equal');
xlabel(ax_main, 'X Position (m)');
ylabel(ax_main, 'Y Position (m)');
title(ax_main, 'Interactive 2D Rocket Trajectory');

%% Plot complete trajectory path (faded)
x_traj = X_data(1, :);
y_traj = X_data(2, :);
plot(ax_main, x_traj, y_traj, 'k--', 'LineWidth', 1);

% Mark landing target
plot(ax_main, P.x_target, P.y_target, 'ro', 'MarkerSize', 12, 'LineWidth', 3);
text(ax_main, P.x_target + 10, P.y_target + 10, 'Target', 'FontSize', 12);

%% Initialize animated elements
% Rocket position marker
h_rocket = plot(ax_main, x_traj(1), y_traj(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% Vector arrows (unit length initially)
vector_scale = 50; % Scaling factor for vector display
h_velocity    = quiver(ax_main, x_traj(1), y_traj(1), 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
h_orientation = quiver(ax_main, x_traj(1), y_traj(1), 0, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5);
h_thrust      = quiver(ax_main, x_traj(1), y_traj(1), 0, 0, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);

% Legend for vectors
legend(ax_main, [h_velocity, h_orientation, h_thrust], ...
       {'Velocity', 'Rocket Up', 'Thrust'}, 'Location', 'northwest');

%% Create control panel
panel_height = 0.12;
panel_y      = 0.02;
button_width = 0.08;
button_height = 0.06;
slider_width = 0.4;

% Control panel background
uipanel('Parent', fig_handle, 'Position', [0.1, panel_y, 0.8, panel_height], ...
        'BackgroundColor', [0.9, 0.9, 0.9]);

% Play/Pause button
btn_play = uicontrol('Parent', fig_handle, 'Style', 'pushbutton', ...
                     'String', 'Play', 'Units', 'normalized', ...
                     'Position', [0.3, panel_y + 0.03, button_width, button_height], ...
                     'FontSize', 12, 'Callback', @play_pause_callback);

% Restart button
btn_restart = uicontrol('Parent', fig_handle, 'Style', 'pushbutton', ...
                        'String', 'Restart', 'Units', 'normalized', ...
                        'Position', [0.4, panel_y + 0.03, button_width, button_height], ...
                        'FontSize', 12, 'Callback', @restart_callback);

% Backward button
btn_backward = uicontrol('Parent', fig_handle, 'Style', 'pushbutton', ...
                         'String', '<<', 'Units', 'normalized', ...
                         'Position', [0.2, panel_y + 0.03, button_width/2, button_height], ...
                         'FontSize', 12, 'Callback', @backward_callback);

% Forward button  
btn_forward = uicontrol('Parent', fig_handle, 'Style', 'pushbutton', ...
                        'String', '>>', 'Units', 'normalized', ...
                        'Position', [0.5, panel_y + 0.03, button_width/2, button_height], ...
                        'FontSize', 12, 'Callback', @forward_callback);

% Timeline slider
slider_timeline = uicontrol('Parent', fig_handle, 'Style', 'slider', ...
                            'Units', 'normalized', ...
                            'Position', [0.55, panel_y + 0.05, slider_width, 0.02], ...
                            'Min', 1, 'Max', n_points, 'Value', 1, ...
                            'SliderStep', [1/(n_points-1), 10/(n_points-1)], ...
                            'Callback', @slider_callback);

% Time display
txt_time = uicontrol('Parent', fig_handle, 'Style', 'text', ...
                     'String', sprintf('t = %.2f s', t_data(1)), ...
                     'Units', 'normalized', ...
                     'Position', [0.55, panel_y + 0.08, 0.15, 0.03], ...
                     'FontSize', 10, 'BackgroundColor', [0.9, 0.9, 0.9]);

%% Animation control variables
current_idx   = 1;
is_playing    = false;
play_timer    = [];
play_speed    = 0.05; % seconds between frames
play_direction = 1;   % 1 for forward, -1 for backward

%% Initialize display
update_display(current_idx);

%% Callback functions
function play_pause_callback(~, ~)
    if is_playing
        % Pause
        stop_animation();
        set(btn_play, 'String', 'Play');
    else
        % Play
        start_animation();
        set(btn_play, 'String', 'Pause');
    end
end

function restart_callback(~, ~)
    stop_animation();
    current_idx = 1;
    set(slider_timeline, 'Value', current_idx);
    update_display(current_idx);
    set(btn_play, 'String', 'Play');
end

function backward_callback(~, ~)
    stop_animation();
    play_direction = -1;
    start_animation();
    set(btn_play, 'String', 'Pause');
end

function forward_callback(~, ~)
    stop_animation();
    play_direction = 1;
    start_animation();
    set(btn_play, 'String', 'Pause');
end

function slider_callback(src, ~)
    stop_animation();
    current_idx = round(get(src, 'Value'));
    update_display(current_idx);
    set(btn_play, 'String', 'Play');
end

function start_animation()
    if ~isempty(play_timer)
        stop(play_timer);
        delete(play_timer);
    end
    
    play_timer = timer('ExecutionMode', 'fixedRate', ...
                       'Period', play_speed, ...
                       'TimerFcn', @animate_step);
    start(play_timer);
    is_playing = true;
end

function stop_animation()
    if ~isempty(play_timer) && isvalid(play_timer)
        stop(play_timer);
        delete(play_timer);
        play_timer = [];
    end
    is_playing = false;
end

function animate_step(~, ~)
    current_idx = current_idx + play_direction;
    
    % Check bounds and stop if reached
    if current_idx > n_points || current_idx < 1
        current_idx = max(1, min(n_points, current_idx));
        stop_animation();
        set(btn_play, 'String', 'Play');
        return;
    end
    
    set(slider_timeline, 'Value', current_idx);
    update_display(current_idx);
end

function update_display(idx)
    % Get current state
    x_pos  = X_data(1, idx);
    y_pos  = X_data(2, idx);
    
    % Update rocket position
    set(h_rocket, 'XData', x_pos, 'YData', y_pos);
    
    % Update vectors
    vel_vec    = velocity_data(:, idx);
    orient_vec = orient_data(:, idx);  
    thrust_vec = thrust_data(:, idx);
    
    % Normalize vectors for display (unit length with scaling)
    vel_norm    = vel_vec / max(norm(vel_vec), 1e-6) * vector_scale;
    orient_norm = orient_vec * vector_scale; % Already unit length
    thrust_norm = thrust_vec / max(norm(thrust_vec), 1e-6) * vector_scale;
    
    % Handle zero vectors gracefully
    if norm(vel_vec) < 1e-6
        vel_norm = [0; 0];
    end
    if norm(thrust_vec) < 1e-6
        thrust_norm = [0; 0];
    end
    
    % Update vector displays
    set(h_velocity, 'XData', x_pos, 'YData', y_pos, 'UData', vel_norm(1), 'VData', vel_norm(2));
    set(h_orientation, 'XData', x_pos, 'YData', y_pos, 'UData', orient_norm(1), 'VData', orient_norm(2));
    set(h_thrust, 'XData', x_pos, 'YData', y_pos, 'UData', thrust_norm(1), 'VData', thrust_norm(2));
    
    % Update camera to follow rocket
    camera_width  = 400; % meters
    camera_height = 300; % meters
    xlim(ax_main, [x_pos - camera_width/2, x_pos + camera_width/2]);
    ylim(ax_main, [y_pos - camera_height/2, y_pos + camera_height/2]);
    
    % Update time display
    set(txt_time, 'String', sprintf('t = %.2f s', t_data(idx)));
    
    % Force redraw
    drawnow;
end

%% Cleanup function when figure is closed
set(fig_handle, 'CloseRequestFcn', @close_callback);

function close_callback(~, ~)
    stop_animation();
    delete(fig_handle);
end

end