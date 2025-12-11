clear; clc; close all;

% Input Function Here
theta = deg2rad(0:720);
r = 1 + 0.5*sin(2*theta);

polarplot(theta, r, 'LineWidth', 1.5);
ax = gca;
ax.ThetaZeroLocation = 'right';
ax.ThetaDir = 'clockwise';
ax.RLim = [0 3];

addCrankLabels(ax, 30, 3.2, 3.5); % tick every 30°, auto radius placement

function addCrankLabels(ax, tickStep, r1, r2)
% addCrankLabels(ax, tickStep, r1, r2)
% Adds 0–360° and 360–720° angle labels to a polar plot without overlap.
%
% ax        – handle to polaraxes (e.g., gca)
% tickStep  – angle step in degrees (e.g., 30 or 45)
% r1, r2    – radii for inner and outer label rings

if nargin < 1 || isempty(ax)
    ax = gca;
end
if nargin < 2, tickStep = 30; end
if nargin < 3, r1 = max(ax.RLim)*1.05; end
if nargin < 4, r2 = max(ax.RLim)*1.15; end

% Disable built-in theta tick labels
ax.ThetaTickLabel = [];

% First revolution (0–360°)
tickAngles = 0:tickStep:330;
for ang = tickAngles
    t = deg2rad(ang);
    text(ax, t, r1, sprintf('%d°', ang),'HorizontalAlignment','center','VerticalAlignment','middle');
end

% Second revolution (360–720°)
for ang = tickAngles
    t = deg2rad(ang);
    text(ax, t, r2, sprintf('%d°', ang+360),'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.4 0.4 0.4]); % lighter gray for clarity
end

% Stack top labels manually (0°, 360°, 720°)
text(ax, 0, r1, '0°', 'HorizontalAlignment','center');
text(ax, 0, r2, '360°', 'HorizontalAlignment','center', 'Color',[0.4 0.4 0.4]);
text(ax, 0, r2 + 0.35, '720°', 'HorizontalAlignment','center');

end