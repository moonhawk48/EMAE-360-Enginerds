% Flat twin, 180° crankpins
% Plots: 
% (1) primary & secondary shaking forces: cyl A, cyl B, and their sum
% (2) rocking couple: total plus primary and secondary components

% Conventions:
% m_rec: reciprocating mass per cylinder (piston + small-end eq. rod)
% r: crank radius; l: conrod length; B: cylinder centerline spacing
% Cylinder A and B produce equal-and-opposite line-of-stroke forces.
% Their vector sum (shaking force) cancels
% Their separation B creates a rocking couple

clear; clc;


% Inputs

r     = 0.0432;            % [m] crank radius (4.32 cm)
l     = 0.1524;          % [m] conrod length (15.24 cm)
m_rec = 2.016141732;     % [kg] reciprocating mass per cylinder
rpm   = 7031;            % [rev/min] evaluation speed
B     = 0.05;            % [m] spacing between cylinder centerlines (measure in CAD)


% Derived and angle grid

omega = rpm * 2*pi/60;             % [rad/s]
n     = l/r;                        % conrod ratio
th    = linspace(0, 2*pi, 2001);    % 0..360 deg inclusive
tdeg  = th * 180/pi;


% One-cylinder components along its line of stroke
% (1× = primary, 2× = secondary)

F1_single = m_rec * r * omega^2 .* cos(th);
F2_single = m_rec * r * omega^2 * (1/n) .* cos(2*th);


% Boxer twin: equal-and-opposite line-of-stroke forces
% Use a common global x-direction
% Cylinder B force is the opposite of A
% This gives zero resultant shaking force.
% The separated lines create a couple
F1_A = +F1_single;
F1_B = -F1_single;
F2_A = +F2_single;
F2_B = -F2_single;

F_shake_primary_sum   = F1_A + F1_B;   % should be 0
F_shake_secondary_sum = F2_A + F2_B;   % should be 0

% Rocking couple: forces act in parallel planes separated by B
% Couple about the crank centerline (z-axis out of the page):
M_primary   = B .* F1_A;              % primary component
M_secondary = B .* F2_A;              % secondary component
M_total     = M_primary + M_secondary;

% For readability
Fscale = 1e-3;   % N -> kN
Mscale = 1e-3;   % N*m -> kN*m


% FIGURE 1: Shaking forces (A, B, and Sum) for primary and secondary

figure('Color','w','Position',[80 80 980 720]);

subplot(2,1,1); % Primary (1x)
plot(tdeg, F1_A*Fscale, 'LineWidth',1.5); hold on;
plot(tdeg, F1_B*Fscale, 'LineWidth',1.5);
plot(tdeg, F_shake_primary_sum*Fscale, 'k','LineWidth',2);
grid on; xlim([0 360]);
xlabel('Crank Angle [deg]'); ylabel('Primary Shaking Force [kN]');
title('Primary (1×) Shaking Force');
legend('Cylinder A','Cylinder B','Sum (A+B)');

subplot(2,1,2); % Secondary (2x)
plot(tdeg, F2_A*Fscale, 'LineWidth',1.5); hold on;
plot(tdeg, F2_B*Fscale, 'LineWidth',1.5);
plot(tdeg, F_shake_secondary_sum*Fscale, 'k','LineWidth',2);
grid on; xlim([0 360]);
xlabel('Crank angle [deg]'); ylabel('Secondary Shaking Force [kN]');
title('Secondary (2×) Shaking Force');
legend('Cylinder A','Cylinder B','Sum (A+B)');

sgtitle(sprintf('Shaking Forces'));

% FIGURE 2: Rocking couple (total as well as components)

figure('Color','w','Position',[100 100 980 420]);
plot(tdeg, M_total*Mscale, 'LineWidth',2); hold on;
plot(tdeg, M_primary*Mscale, 'LineWidth',1.5);
plot(tdeg, M_secondary*Mscale, 'LineWidth',1.5);
grid on; xlim([0 360]);
xlabel('Crank angle [deg]'); ylabel('Rocking Couple [kN·m]');
title(sprintf('Rocking Couple About Crank Centerline at %0.0f rpm', rpm));
legend('Total (primary + secondary)','Primary component','Secondary component','Location','southoutside','Orientation','horizontal');

% Quick numbers
F1_pk = max(abs(F1_single));
F2_pk = max(abs(F2_single));
M1_pk = B * F1_pk;
M2_pk = B * F2_pk;
fprintf('Primary shaking sum (should be ~0): %g N peak\n', max(abs(F_shake_primary_sum)));
fprintf('Secondary shaking sum (should be ~0): %g N peak\n', max(abs(F_shake_secondary_sum)));
fprintf('Primary rocking couple peak: %g N·m\n', M1_pk);
fprintf('Secondary rocking couple peak: %g N·m\n', M2_pk);
