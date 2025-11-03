clear; clc; close all;

% Input parameters
Lv_max = 16.9;              % Max valve lift [mm]
beta = 90;                  % Lift angle [deg]
rpm = 2500;                 % Camshaft speed [rev/min]
omega = rpm * 2*pi/60;      % Angular speed [rad/s]

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
theta = u * beta;           % Crank angle [deg]

% Base equations
S = a4*u.^4 + a5*u.^5 + a6*u.^6 + a7*u.^7;                        % Displacement [mm]
dS_du = 4*a4*u.^3 + 5*a5*u.^4 + 6*a6*u.^5 + 7*a7*u.^6;            % First derivative
d2S_du2 = 12*a4*u.^2 + 20*a5*u.^3 + 30*a6*u.^4 + 42*a7*u.^5;      % Second derivative
d3S_du3 = 24*a4*u + 60*a5*u.^2 + 120*a6*u.^3 + 210*a7*u.^4;       % Third derivative

% Convert to real units (mm/s, mm/s², mm/s³)
V = dS_du * (omega/beta*pi/180);             % Velocity [mm/s]
A = d2S_du2 * (omega/beta*pi/180)^2;         % Acceleration [mm/s²]
J = d3S_du3 * (omega/beta*pi/180)^3;         % Jerk [mm/s³]

% --- Scale x-axis to degrees ---
theta1 = u * 90;  % 0–90 degrees

% --- Mirror motion to 90–180° ---
theta2 = 180 - fliplr(theta1);  % 90–180 degrees
Lv2 = fliplr(S);
Lv_dot2 = -fliplr(V);     % velocity reverses
Lv_ddot2 = fliplr(A);     % acceleration same direction
Lv_tdot2 = -fliplr(J);    % jerk reverses

% --- Combine lift + return ---
theta = [theta1, theta2];
S = [S, Lv2];
V = [V, Lv_dot2];
A = [A, Lv_ddot2];
J = [J, Lv_tdot2];


% --- Plot SVAJ ---
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