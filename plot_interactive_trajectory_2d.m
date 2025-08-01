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
clf(fig_handle);

%% Setup main axes
ax_main = axes('Parent', fig_handle, 'Position', [0.1, 0.25, 0.8, 0.65]);
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
% Load rocket image
rocket_image_path = './figures_2d/f9_image.png';
try
    [rocket_img, ~, alpha] = imread(rocket_image_path);
    has_rocket_image = true;
catch
    warning('Could not load rocket image from %s, using marker instead', rocket_image_path);
    has_rocket_image = false;
end

% Rocket position marker (red dot for center of gravity)
h_rocket = plot(ax_main, x_traj(1), y_traj(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% Rocket image will be created in update_display function

% Vector arrows (unit length initially)
vector_scale = 50; % Scaling factor for vector display
h_velocity    = quiver(ax_main, x_traj(1), y_traj(1), 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
h_orientation = quiver(ax_main, x_traj(1), y_traj(1), 0, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5);
h_thrust      = quiver(ax_main, x_traj(1), y_traj(1), 0, 0, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);

% Legend for vectors
legend(ax_main, [h_velocity, h_orientation, h_thrust], ...
       {'Velocity', 'Rocket Up', 'Thrust'}, 'Location', 'northwest');

%% Create control panel
panel_height = 0.18;
panel_y      = 0.02;
button_width = 0.08;
button_height = 0.05;
slider_width = 0.35;

% Control panel background
uipanel('Parent', fig_handle, 'Position', [0.05, panel_y, 0.9, panel_height], ...
        'BackgroundColor', [0.9, 0.9, 0.9]);

% Center the buttons better
btn_center_x = 0.25 - (2*button_width + button_width)/2; % Center point for 4 buttons

% Backward button
btn_backward = uicontrol('Parent', fig_handle, 'Style', 'pushbutton', ...
                         'String', '<<', 'Units', 'normalized', ...
                         'Position', [btn_center_x - 0.05, panel_y + 0.09, button_width/2, button_height], ...
                         'FontSize', 12, 'Callback', @backward_callback);

% Play/Pause button
btn_play = uicontrol('Parent', fig_handle, 'Style', 'pushbutton', ...
                     'String', 'Play', 'Units', 'normalized', ...
                     'Position', [btn_center_x + 0.01, panel_y + 0.09, button_width, button_height], ...
                     'FontSize', 12, 'Callback', @play_pause_callback);

% Restart button
btn_restart = uicontrol('Parent', fig_handle, 'Style', 'pushbutton', ...
                        'String', 'Restart', 'Units', 'normalized', ...
                        'Position', [btn_center_x + 0.10, panel_y + 0.09, button_width, button_height], ...
                        'FontSize', 12, 'Callback', @restart_callback);

% Forward button  
btn_forward = uicontrol('Parent', fig_handle, 'Style', 'pushbutton', ...
                        'String', '>>', 'Units', 'normalized', ...
                        'Position', [btn_center_x + 0.19, panel_y + 0.09, button_width/2, button_height], ...
                        'FontSize', 12, 'Callback', @forward_callback);

% Timeline slider (thicker)
slider_timeline = uicontrol('Parent', fig_handle, 'Style', 'slider', ...
                            'Units', 'normalized', ...
                            'Position', [0.46, panel_y + 0.11, slider_width, 0.025], ...
                            'Min', 1, 'Max', n_points, 'Value', 1, ...
                            'SliderStep', [1/100, 10/100], ...
                            'Callback', @slider_callback);
%                             'SliderStep', [1/(n_points-1), 10/(n_points-1)], ...

% Timeline labels
txt_time_start = uicontrol('Parent', fig_handle, 'Style', 'text', ...
                          'String', sprintf('%.1fs', t_data(1)), ...
                          'Units', 'normalized', ...
                          'Position', [0.4, panel_y + 0.105, 0.05, 0.025], ...
                          'FontSize', 9, 'BackgroundColor', [0.9, 0.9, 0.9], ...
                          'HorizontalAlignment', 'right');

txt_time_end = uicontrol('Parent', fig_handle, 'Style', 'text', ...
                        'String', sprintf('%.1fs', t_data(end)), ...
                        'Units', 'normalized', ...
                        'Position', [0.82, panel_y + 0.105, 0.05, 0.025], ...
                        'FontSize', 9, 'BackgroundColor', [0.9, 0.9, 0.9], ...
                        'HorizontalAlignment', 'left');

% Current time display
txt_time = uicontrol('Parent', fig_handle, 'Style', 'text', ...
                     'String', sprintf('t = %.2f s', t_data(1)), ...
                     'Units', 'normalized', ...
                     'Position', [0.56, panel_y + 0.14, 0.15, 0.03], ...
                     'FontSize', 11, 'BackgroundColor', [0.9, 0.9, 0.9], ...
                     'FontWeight', 'bold');

% Speed control slider (log scale)
speed_min_log = log10(0.1);  % 0.1x speed
speed_max_log = log10(10);   % 10x speed
speed_default_log = log10(1); % 1x speed

slider_speed = uicontrol('Parent', fig_handle, 'Style', 'slider', ...
                        'Units', 'normalized', ...
                        'Position', [0.46, panel_y + 0.04, slider_width, 0.025], ...
                        'Min', speed_min_log, 'Max', speed_max_log, 'Value', speed_default_log, ...
                        'SliderStep', [0.05, 0.2], ...
                        'Callback', @speed_slider_callback);

% Speed control label
% txt_speed_label = uicontrol('Parent', fig_handle, 'Style', 'text', ...
%                            'String', 'Speed:', ...
%                            'Units', 'normalized', ...
%                            'Position', [0.5, panel_y + 0.05, 0.06, 0.025], ...
%                            'FontSize', 10, 'BackgroundColor', [0.9, 0.9, 0.9], ...
%                            'HorizontalAlignment', 'right');

% Speed graduations
speed_values = [0.1, 0.2, 0.5, 1, 2, 5, 10];
for i = 1:length(speed_values)
    speed_val = speed_values(i);
    slider_pos = (log10(speed_val) - speed_min_log) / (speed_max_log - speed_min_log);
    x_pos = 0.46 + slider_pos * slider_width;
    
    % Graduation mark
    uicontrol('Parent', fig_handle, 'Style', 'text', ...
              'String', '|', ...
              'Units', 'normalized', ...
              'Position', [x_pos - 0.005, panel_y + 0.065, 0.01, 0.015], ...
              'FontSize', 8, 'BackgroundColor', [0.9, 0.9, 0.9], ...
              'HorizontalAlignment', 'center');
    
    % Graduation label
    if speed_val < 1
        label_str = sprintf('%.1fx', speed_val);
    else
        label_str = sprintf('%.0fx', speed_val);
    end
    uicontrol('Parent', fig_handle, 'Style', 'text', ...
              'String', label_str, ...
              'Units', 'normalized', ...
              'Position', [x_pos - 0.015, panel_y + 0.08, 0.03, 0.015], ...
              'FontSize', 8, 'BackgroundColor', [0.9, 0.9, 0.9], ...
              'HorizontalAlignment', 'center');
end

% Current speed display
txt_speed = uicontrol('Parent', fig_handle, 'Style', 'text', ...
                     'String', 'Speed: 1.0x', ...
                     'Units', 'normalized', ...
                     'Position', [0.58, panel_y + 0.01, 0.1, 0.025], ...
                     'FontSize', 10, 'BackgroundColor', [0.9, 0.9, 0.9], ...
                     'FontWeight', 'bold');

%% Animation control variables
current_idx      = 1.0;  % Now supports fractional indices
is_playing       = false;
play_timer       = [];
speed_multiplier = 1.0;  % Current speed multiplier (1.0 = real time)
display_fps      = 30;   % Fixed display frame rate
play_direction   = 1;    % 1 for forward, -1 for backward

% Rocket image handle (declared here for proper scope)
h_rocket_image   = [];

% Fixed timer period and simulation parameters
timer_period = 1 / display_fps;  % Fixed 33ms timer period
dt_sim_avg = mean(diff(t_data)); % Average time step in simulation data
steps_per_frame = (timer_period * speed_multiplier) / dt_sim_avg; % Steps to advance per timer tick

%% Initialize display
update_display(current_idx);

pause(0.5)
plotedit(fig_handle,'off')


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
    current_idx = 1.0;
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
    current_idx = get(src, 'Value'); % No rounding - keep fractional values
    update_display(current_idx);
    set(btn_play, 'String', 'Play');
end

function speed_slider_callback(src, ~)
    % Get log value and convert to linear speed multiplier
    log_val = get(src, 'Value');
    speed_multiplier = 10^log_val;
    
    % Recalculate steps per frame with new speed
    steps_per_frame = (timer_period * speed_multiplier) / dt_sim_avg;
    
    % Update speed display
    set(txt_speed, 'String', sprintf('Speed: %.1fx', speed_multiplier));
    
    % Timer period remains fixed - no need to restart timer
    % The animate_step function will use the updated steps_per_frame automatically
end

function start_animation()
    if ~isempty(play_timer)
        stop(play_timer);
        delete(play_timer);
    end
    
    play_timer = timer('ExecutionMode', 'fixedRate', ...
                       'Period', timer_period, ...
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
    % Calculate step size based on current speed and direction
    step_size = steps_per_frame * play_direction;
    current_idx = current_idx + step_size;
    
    % Check bounds and stop if reached
    if current_idx > n_points || current_idx < 1.0
        current_idx = max(1.0, min(n_points, current_idx));
        stop_animation();
        set(btn_play, 'String', 'Play');
        return;
    end
    
    update_display(current_idx);
end

function [interp_state, interp_vectors, interp_time] = interpolate_data(fractional_idx)
    % Interpolate between simulation data points for fractional indices
    
    % Clamp to valid range
    fractional_idx = max(1.0, min(n_points, fractional_idx));
    
    % Get integer indices for interpolation
    idx_low = floor(fractional_idx);
    idx_high = min(idx_low + 1, n_points);
    
    % Calculate interpolation weight
    if idx_low == idx_high
        weight = 0; % Exact integer index
    else
        weight = fractional_idx - idx_low;
    end
    
    % Interpolate state data
    state_low = X_data(:, idx_low);
    state_high = X_data(:, idx_high);
    interp_state = state_low + weight * (state_high - state_low);
    
    % Interpolate vector data
    vel_low = velocity_data(:, idx_low);
    vel_high = velocity_data(:, idx_high);
    orient_low = orient_data(:, idx_low);
    orient_high = orient_data(:, idx_high);
    thrust_low = thrust_data(:, idx_low);
    thrust_high = thrust_data(:, idx_high);
    
    interp_vectors.velocity = vel_low + weight * (vel_high - vel_low);
    interp_vectors.orientation = orient_low + weight * (orient_high - orient_low);
    interp_vectors.thrust = thrust_low + weight * (thrust_high - thrust_low);
    
    % Interpolate time
    time_low = t_data(idx_low);
    time_high = t_data(idx_high);
    interp_time = time_low + weight * (time_high - time_low);
end

function update_display(fractional_idx)
    % Get interpolated data for current fractional index
    [interp_state, interp_vectors, interp_time] = interpolate_data(fractional_idx);
    
    % Extract position and orientation
    x_pos = interp_state(1);
    y_pos = interp_state(2);
    theta = interp_state(5); % Rocket pitch angle
    
    % Update rocket position (center of gravity marker)
    set(h_rocket, 'XData', x_pos, 'YData', y_pos);
    
    % Update rocket image if available
    if has_rocket_image
        % Delete previous image if it exists
        if ~isempty(h_rocket_image) && isvalid(h_rocket_image)
            delete(h_rocket_image);
        end
        
        % Calculate image size and position
        image_scale = 0.05;
        img_height = size(rocket_img, 1) * image_scale;
        img_width = size(rocket_img, 2) * image_scale;
        
        % Create rotated image
        % Rotate image by theta (rocket pitch angle)
        % Note: imrotate rotates counterclockwise, theta is measured from vertical
        rotated_img = imrotate(rocket_img, -rad2deg(theta) + 180, 'bilinear', 'crop');
        rotated_alpha = imrotate(alpha, -rad2deg(theta) + 180, 'bilinear', 'crop');
        
        % Calculate image extents centered on rocket position
        x_extent = [x_pos - img_width/2, x_pos + img_width/2];
        y_extent = [y_pos - img_height/2, y_pos + img_height/2];
        
        % Display the rotated image
        h_rocket_image = image(ax_main, x_extent, y_extent, rotated_img);
        set(h_rocket_image, 'AlphaData', rotated_alpha);
        
        % Ensure rocket marker and vectors are on top of the image
        uistack(h_rocket, 'top');
        uistack(h_velocity, 'top');
        uistack(h_orientation, 'top');
        uistack(h_thrust, 'top');
    end
    
    % Get interpolated vectors
    vel_vec    = interp_vectors.velocity;
    orient_vec = interp_vectors.orientation;  
    thrust_vec = interp_vectors.thrust;
    
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
    set(txt_time, 'String', sprintf('t = %.2f s', interp_time));
    
    % Update timeline slider position
    set(slider_timeline, 'Value', fractional_idx);
    
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