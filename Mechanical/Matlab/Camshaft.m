clear; clc; close all;

% Input parameters
Lv_max = 16.9; % Max valve lift [mm]
beta = 90; % Lift angle [deg]
rpm = 2500; % Camshaft speed [rev/min]
omega = rpm * 2*pi/60; % Angular speed [rad/s]

% Solving for 4-5-6-7 Polynomial Coefficients
b = 1;
A = [b^4 b^5 b^6 b^7; 4*b^3 5*b^4 6*b^5 7*b^6; 12*b^2 20*b^3 30*b^4 42*b^5; 24*b 60*b^2 120*b^3 210*b^4];
B = [Lv_max; 0; 0; 0];
coefficients = A \ B;
a4 = coefficients(1);
a5 = coefficients(2);
a6 = coefficients(3);
a7 = coefficients(4);

% Normalized motion variable
u = 0:0.001:1;
theta = u * beta; % Crank angle [deg]

% Base equations
S = a4*u.^4 + a5*u.^5 + a6*u.^6 + a7*u.^7; % Displacement [mm]
dS_du = 4*a4*u.^3 + 5*a5*u.^4 + 6*a6*u.^5 + 7*a7*u.^6; % Velocity
d2S_du2 = 12*a4*u.^2 + 20*a5*u.^3 + 30*a6*u.^4 + 42*a7*u.^5; % Acceleration
d3S_du3 = 24*a4*u + 60*a5*u.^2 + 120*a6*u.^3 + 210*a7*u.^4; % Jerk

% Convert to units (mm/s, mm/s², mm/s³)
V = dS_du * (omega/beta*pi/180); % Velocity [mm/s]
A = d2S_du2 * (omega/beta*pi/180)^2; % Acceleration [mm/s²]
J = d3S_du3 * (omega/beta*pi/180)^3; % Jerk [mm/s³]

% Scale x-axis to degrees
theta1 = u * 90;  % 0–90 degrees

% Mirror motion from 90–180°
theta2 = 180 - fliplr(theta1);  % 90–180 degrees
Lv2 = fliplr(S);
Lv_dot2 = -fliplr(V); % velocity reverses
Lv_ddot2 = fliplr(A); % acceleration remains the same
Lv_tdot2 = -fliplr(J); % jerk reverses

% Combining lift + return
theta = [theta1, theta2];
S = [S, Lv2];
V = [V, Lv_dot2];
A = [A, Lv_ddot2];
J = [J, Lv_tdot2];

% Plot SVAJ
figure;
sgtitle(sprintf('SVAJ Curves Using 4-5-6-7 Polynomials @ %d RPM', rpm));

subplot(4,1,1);
plot(theta, S, 'LineWidth', 1.5);
xlabel('Crank Angle θ [deg]');
ylabel('S [mm]');
grid on;

subplot(4,1,2);
plot(theta, V, 'LineWidth', 1.5);
xlabel('Crank Angle θ [deg]');
ylabel('V [mm/s]');
grid on;

subplot(4,1,3);
plot(theta, A, 'LineWidth', 1.5);
xlabel('Crank Angle θ [deg]');
ylabel('A [mm/s²]');
grid on;

subplot(4,1,4);
plot(theta, J, 'LineWidth', 1.5);
xlabel('Crank Angle θ [deg]');
ylabel('J [mm/s³]');
grid on;

% Cam Profile for Flat-Faced Follower
R_b = 10; % Base circle radius [mm]

% Convert S and dS/dθ to proper units
theta_total = linspace(0, 360, length(S));  % 0–360 degrees SEPARATE INTO 90-90-180
theta_rad = deg2rad(theta_total);
dS_dtheta = [diff(S)./diff(theta_rad), 0];  % mm/rad

% Flat-faced follower contact coordinates
x = (R_b + S) .* cos(theta_rad) - dS_dtheta .* sin(theta_rad);
y = (R_b + S) .* sin(theta_rad) + dS_dtheta .* cos(theta_rad);

% Plot cam profile
figure;
plot(x, y, 'r', 'LineWidth', 1.5); hold on;

% Base circle for reference
ang = linspace(0, 2*pi, 200);
plot(R_b*sin(ang), R_b*cos(ang), 'b--', 'LineWidth', 1.5);

% Formatting
axis equal;
xlim([-30,30]);
ylim([-30,30]);
xlabel('X [mm]');
ylabel('Y [mm]');
title('Cam Profile');
grid on;

% Exporting Profile to Solidworks
x_col = x(:);
y_col = y(:);
z_col = zeros(size(x_col));  % flat in XY plane

% Combine into one matrix
cam_profile_xyz = [x_col, y_col, z_col];

% Save as CSV
writematrix(cam_profile_xyz, 'cam_profile_points.csv');
%%

% Cam Profile Mapping
R_b = 5; % Base circle radius [mm]

% Add base radius and extend motion to full 360°
S_total = [S, zeros(1, 180)];  % add dwell back to base circle
theta_total = [linspace(0, 180, length(S)), linspace(180, 360, 180)];
theta_rad = deg2rad(theta_total);

% Cam polar coordinates
r = R_b + S_total;
x = r .* cos(theta_rad);
y = r .* sin(theta_rad);

% Plot Profile
figure(2);
plot(x, y, 'LineWidth', 1.5);
hold on;

% Plot Base Circle
angle = linspace(0, 2*pi, 200);
x_circle = R_b * cos(angle);
y_circle = R_b * sin(angle);
plot(x_circle, y_circle, 'b--', 'LineWidth', 1.5);

% Center and formatting
plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
xlim([-30, 30]);
ylim([-30, 30]);
title('Cam Lobe Profile');
xlabel('X [mm]');
ylabel('Y [mm]');
grid on;

%%
% Cam Profile Mapping
R_b = 5; % Base circle radius [mm]
S = [S + R_b, R_b*ones(1,length(theta))]; % Adding Base Radius to Lift from 180 to 360 deg.
V = [V, zeros(1,length(theta))]; % Adding zero velocity from 180 to 360 deg.
u = 1:0.0005:2;
theta3 = u*2*beta;
theta = deg2rad([theta, theta3]); % Crank Angle [rad]

% Mapping Profile
for i = 1:1:length(theta)
    x(i) = (R_b + S(i))*cos(theta(i)) - V(i)*sin(theta(i)); % X Function
    y(i) = (R_b + S(i))*sin(theta(i)) + V(i)*cos(theta(i)); % Y Function
end

% Plot Profile
figure(2);
plot(x,y,'LineWidth',1.5);

% Plot Base Circle
hold on
radius = R_b;
xCenter = 0;
yCenter = 0;
angle = linspace(0, 2*pi, 100); % Angles from 0 to 2*pi, with 100 points
x_circle = radius * cos(angle) + xCenter;
y_circle = radius * sin(angle) + yCenter;
plot(x_circle, y_circle, 'b--', 'LineWidth', 2); % 'b-' for a blue solid line
hold on
% Optionally, plot the center of the circle
plot(xCenter, yCenter, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
xlim([-30,30]);
ylim([-30,30]);