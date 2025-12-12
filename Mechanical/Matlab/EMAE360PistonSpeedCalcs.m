% PISTON SPEED PLOTS people e
% Uses exact slider–crank position and velocity formulas.
% Plots:
% Max piston speed vs RPM over a chosen RPM vector
% Piston speed vs crank angle at a chosen RPM

clear; 
clc;

% Geometry (Need to adjust these as CAD is finalized)

r = 0.0432; % [m] crank radius  (5 cm)
l = 0.18; % [m] conrod length (15.24 cm)

% RPM settings

rpm_vec = 800:200:6000; % range for max-speed curve (may need adjustments)
rpm_plot = 5000; % single RPM for speed-vs-angle plot

% Angle from 0 to 360 degrees

th = linspace(0, 2*pi, 4001);
tdeg = th * 180/pi;

% Exact slider–crank kinematics
% Piston position from crank centerline along line of stroke:
% x(θ) = r*cosθ + sqrt(l^2 - (r*sinθ)^2)
% Velocity: v(θ) = ω * dx/dθ, with
% dx/dθ = -r*sinθ - (r^2 * sinθ * cosθ) / sqrt(l^2 - r^2*sin^2θ)

dx_dtheta = @(theta) (-r.*sin(theta)) - ( (r.^2 .* sin(theta) .* cos(theta)) ./ sqrt(l.^2 - (r.^2 .* sin(theta).^2) ) );

% Max piston speed vs RPM

vmax = zeros(size(rpm_vec));
for k = 1:numel(rpm_vec)
    omega = rpm_vec(k) * 2*pi/60;
    v = omega .* dx_dtheta(th);    % [m/s]
    vmax(k) = max(abs(v));
end

figure('Color','w','Position',[100 100 900 380]);
plot(rpm_vec, vmax, 'LineWidth',1.8);
grid on;
xlabel('Engine speed [rpm]');
ylabel('Max piston speed |v|_{max} [m/s]');
title(sprintf('Max piston speed vs RPM '));


% Piston speed vs crank angle at a chosen RPM

omega_plot = rpm_plot * 2*pi/60;
v_plot = omega_plot .* dx_dtheta(th); % [m/s]

figure('Color','w','Position',[100 520 900 380]);
plot(tdeg, v_plot, 'LineWidth',1.8);
grid on; xlim([0 360]);
xlabel('Crank angle [deg]');
ylabel('Piston speed v(θ) [m/s]');
title(sprintf('Piston speed vs crank angle at %d rpm', rpm_plot));

% Quick readouts
fprintf('Rod ratio l/r = %.3f\n', l/r);
fprintf('At %d rpm: max |v| = %.3f m/s\n', rpm_plot, max(abs(v_plot)));
