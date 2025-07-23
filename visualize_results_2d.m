function visualize_results_2d(hist, P, fm)
%VISUALIZE_RESULTS_2D Plots the results of the 2D MPC simulation.

t = hist.t;
X = hist.X;
U = hist.U;

% Unpack states for clarity
% State unpacking: [x; y; vx; vy; theta; omega; m]
x = X(1,:); y = X(2,:);
vx = X(3,:); vy = X(4,:);
th = X(5,:); om = X(6,:);
% Compute angle of attack (alpha) at each timestep
alpha = zeros(size(t));
for i = 1:length(t)
    vx_i = vx(i); vy_i = vy(i); th_i = th(i);
    V_sq = vx_i^2 + vy_i^2;
    if V_sq > 1e-6
        v_I = [vx_i; vy_i];
        v_aero_I = -v_I / sqrt(V_sq);
        T_B_I = [cos(th_i), sin(th_i); -sin(th_i), cos(th_i)];
        v_aero_B = T_B_I * v_aero_I;
        alpha(i) = atan2(v_aero_B(1), v_aero_B(2));
    else
        alpha(i) = 0;
    end
end
m = X(7,:);
T = U(1,:); delta = U(2,:);

% --- Main 2D Trajectory Plot ---
fm.newFigure('MPC_2D_Trajectory');
plot(x, y, 'b-', 'LineWidth', 2);
hold on;
% Plot rocket orientation at intervals
indices = round(linspace(1, length(t), 15));
for idx = indices
    L = P.L_rocket / 4; % Scaled length for plotting
    px = x(idx); py = y(idx); 
    theta = th(idx);
    
    % Rocket body frame vectors in inertial frame
    % T_I_B = [cos(theta), -sin(theta); sin(theta), cos(theta)]
    % zb-axis (nose direction): T_I_B * [0; 1]
    % xb-axis (side direction): T_I_B * [1; 0]
    T_I_B = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    nose_vec = T_I_B * [0; 1];  % zb unit vector in inertial frame (rocket longitudinal axis)
    side_vec = T_I_B * [1; 0];  % xb unit vector in inertial frame
    
    % Draw rocket body (along zb-axis from tail to nose)
    nose_x = px + L * nose_vec(1);
    nose_y = py + L * nose_vec(2);
    tail_x = px - L * nose_vec(1);
    tail_y = py - L * nose_vec(2);
    plot([tail_x, nose_x], [tail_y, nose_y], 'r-', 'LineWidth', 1.5);
    
    % Draw small cross at nose to show orientation
    cross_size = L/4;
    plot([nose_x - cross_size*side_vec(1), nose_x + cross_size*side_vec(1)], ...
         [nose_y - cross_size*side_vec(2), nose_y + cross_size*side_vec(2)], 'k-', 'LineWidth', 1);
end
plot(P.x0, P.y0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(P.x_target, P.y_target, 'rx', 'MarkerSize', 15, 'LineWidth', 3);
grid on; axis equal;
xlabel('Downrange (m)'); ylabel('Altitude (m)');
title('2D Landing Trajectory with Rocket Attitude');
legend('Trajectory', 'Rocket Attitude', 'Start', 'Target');

% --- Multi-Panel Summary Plot ---
fh = fm.newFigure('MPC_2D_Summary');
fig = fh.Parent;
TL = tiledlayout(fig, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(TL, '2D MPC Simulation Results');

% Position
nexttile;
plot(t, x, 'LineWidth', 1.5); grid on; ylabel('x (m)');
nexttile;
plot(t, y, 'LineWidth', 1.5); grid on; ylabel('y (m)');

% Velocity
nexttile;
plot(t, vx, 'LineWidth', 1.5); grid on; ylabel('v_x (m/s)');
nexttile;
plot(t, vy, 'LineWidth', 1.5); grid on; ylabel('v_y (m/s)');

% Attitude
nexttile;
plot(t, rad2deg(th), 'LineWidth', 1.5); grid on; ylabel('\theta (deg)');
nexttile;
plot(t, rad2deg(alpha), 'LineWidth', 1.5); grid on; ylabel('\alpha (deg)');

% Controls & Mass
nexttile;
plot(t, T/1e3, 'LineWidth', 1.5); grid on;
hold on; yline(P.T_min/1e3, 'k--'); yline(P.T_max/1e3, 'k--');
xlabel('Time (s)'); ylabel('Thrust (kN)');
nexttile;
plot(t, rad2deg(delta), 'LineWidth', 1.5); grid on;
hold on; yline(rad2deg(P.delta_max), 'k--'); yline(rad2deg(-P.delta_max), 'k--');
xlabel('Time (s)'); ylabel('Gimbal \delta (deg)');

% --- Additional Figure: Speed Norm, Mass, Theta, Alpha ---
fh2 = fm.newFigure('MPC_2D_States_Overview');
fig2 = fh2.Parent;
TL2 = tiledlayout(fig2, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(TL2, 'Key States Overview');

% 1. Norm of the speed
nexttile;
plot(t, sqrt(vx.^2 + vy.^2), 'LineWidth', 1.5); grid on;
ylabel('Speed Norm (m/s)');

% 2. Mass
nexttile;
plot(t, m, 'LineWidth', 1.5); grid on;
ylabel('Mass (kg)');

% 3. Theta
nexttile;
plot(t, rad2deg(th), 'LineWidth', 1.5); grid on;
ylabel('\theta (deg)');

% 4. Alpha
nexttile;
plot(t, rad2deg(alpha), 'LineWidth', 1.5); grid on;
ylabel('\alpha (deg)');
xlabel('Time (s)');

end