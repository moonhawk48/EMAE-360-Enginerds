close all
clear

bore = 108; % mm
stroke = 86.4; % mm
l = 180; % mm conrod length
m_cs = 5.62343; % kg mass of crankshaft web (w/ counterweights) for one piston
m_cr = 1.2024; % kg mass of connecting rod
m_p = 0.952; % kg mass of piston head (including wrist pin)
rpm_min = 800;
rpm_max = 6000;
int_rpm = 50;

r = stroke/2; % mm crank radius
r_cs = -9.09; % mm radial distance from rotational axis to crankshaft center of mass
r_cr = 116.69; % mm distance from wrist pin to connecting rod center of mass

m_cp = m_cs*(r_cs/r) + m_cr*(r_cr/l); % approximate mass acting on crankpin
m_wp = m_p + m_cr*((l-r_cr)/l); % approximate mass acting on wrist pin

I = [0.0002 0.0008 0.0203 0.0203 1.83]; % kg m^2 moments of inertia for timing connection, crankpin 1, crankpin 2, and flywheel
% timing inertia is derived from Eq. 3 in the paper, summing all four camshaft moments of inertia at half the crankshaft speed
kt = [1100000 1218040 1218040 1508040]; % Nm/rad from solidworks evaluation
ca = [0 0 1 1 0]; % Nm/rad where pos 1 and 4 should remain 0, 2 and 3 depend on crankpin evaluation
cr = [0 0 0 0]; % Nm/rad relative damping coefficients depend on evaluation of crankshaft

% setup pressure and rpm matrices
opts = detectImportOptions('CylPressurePerRPM.csv'); % pressure in kPa from ?? angle
Raw = readmatrix('CylPressurePerRPM.csv',opts);

rpm_input = Raw(1,2:end);
inc = 720/(length(Raw)-2); % degree increment of data, assuming it includes duplicate values at ends and RPM label row
P = [Raw(2+(170+375)/inc:end,2:end) ; Raw(3:2+(170+375)/inc,2:end)]/1000; % pressures in MPa for crank angles
% angle appears to be offset from spark by 180 degrees
alpha(:,1) = (0:inc:720)*pi/180; % crank angle in radians

% Inertia Matrix
Ialt = m_wp*(r/1000)^2*(1/2 + (r/l)^2/8); % equivallent oscillating inertia
I(3:4) = I(3:4)+Ialt; % adding equivallent inertia to crankpin indecies
M = diag(I);

% Stiffness Matrix
Kt = [kt(1)   -kt(1)         0           0         0 
     -kt(1) kt(1)+kt(2)   -kt(2)         0         0
        0     -kt(2)    kt(2)+kt(3)   -kt(3)       0
        0        0        -kt(3)    kt(3)+kt(4) -kt(4)
        0        0           0        -kt(4)     kt(4)];

% Damping Matrices
Ca = diag(ca);
Cr = [cr(1)   -cr(1)         0           0         0
     -cr(1) cr(1)+cr(2)   -cr(2)         0         0
        0     -cr(2)    cr(2)+cr(3)   -cr(3)       0
        0        0         -cr(3)   cr(3)+cr(4) -cr(4)
        0        0           0        -cr(4)     cr(4)];
C = Cr + Ca;

% First State Matrix
A = [zeros(5,5) eye(5) ; -M\Kt -M\C];


rpm = rpm_min:int_rpm:rpm_max;
theta = zeros(length(rpm),5);
amp = zeros(length(rpm),24,5);
T = zeros(length(rpm),2);
PressureInterpolation;

for rpm_num = 1:length(rpm)
    w = rpm(rpm_num)*pi/30; % rad/s for 5000 rpm
    Fg = (P(:,rpm_num)*pi*bore^2)/4; % N gas load on piston
    beta = asin(sin(alpha)*r/l); % angle between cylider axis and conrod
    Ftg = Fg.*sin(alpha+beta)./cos(beta); % N tangential gas load on crankshaft
    Fi = -m_wp*r*(w^2)*(cos(alpha) + (r/l)*cos(2*alpha) - (r/l)^3*cos(4*alpha)/4 + 9*(r/l)^5*cos(6*alpha)/128)/1000; % N inertial force on wrist pin
    Fti = Fi.*sin(alpha+beta)./cos(beta); % N tangential intertial force on crankshaft
    Ft = Ftg + Fti; % N total tangential load on crankshaft
    Mt = Ft*r/1000; % N.m tangential moment

    %
    % plot Tangential forces vs crankshaft angle for reference
    if (rpm(rpm_num) == 5000)
        figure;
        plot(0:inc:720, Fti/1000, 'b-', 0:inc:720, Ftg/1000, 'r-', 0:inc:720, Ft/1000, 'k--');
        title("Tangential Forces vs Crank Angle at 5000 rpm");
        ylabel("Force (kN)");
        xlabel("Crank Angle (degrees)");
        legend("Inertial Load","Gas Load", "Total Load");
        xticks(0:60:720);
        xlim([0, 720]);
    end
    %}
    %
    % plot instantaneous torque vs crankshaft angle for reference
    if (rpm(rpm_num) == 5000)
        figure;
        plot(0:inc:720, Mt, 'b-', 0:inc:720, [Mt(360/inc:end) ; Mt(2:360/inc)], 'r-'); % second cylinder offset by 360 degrees
        title("Instantaneous Torque vs Crank Angle at 5000 rpm");
        ylabel("Torque (Nm)");
        xlabel("Crank Angle (degrees)");
        legend("Cylinder 1","Cylinder 2");
        xticks(0:60:720);
        xlim([0, 720]);
    end
    %}

    % Fourier expansion of torque
    Cn = 1/length(Mt)*fft(Mt); % fourier coefficients for Mt torque
    % Mt is offset to just TDC before spark for proper values
    i = sqrt(-1); % store imaginary number value for complex number calcs
    Mt(:,1:2) = Cn(1); % start with first coefficient, (average of given Mt)
    for n=1:24 % recursively add subsequent fourier terms
        Mt(:,1) = Mt(:,1) + (Cn(n+1)*exp(i*(n/2)*alpha) + conj(Cn(n+1))*exp(-i*(n/2)*alpha)); % Mt fourier expansion, as per Eq. 17 of paper
        Mt(:,2) = Mt(:,2) + (Cn(n+1)*exp(i*(n/2)*(alpha-2*pi)) + conj(Cn(n+1))*exp(-i*(n/2)*(alpha-2*pi))); % Mt fourier expansion for second piston, offet by 360 degrees
    end
    
    
    %{
    % plot torque fourier expansion vs crankshaft angle for reference
    if (rpm(rpm_num) == 5000)
        figure;
        plot(0:inc:720, Mt(:,1), 'b-', 0:inc:720, Mt(:,2), 'r-');
        title("Fourier Torque vs Crank Angle at 5000 rpm");
        ylabel("Torque (Nm)");
        xlabel("Crank Angle (degrees)");
        legend("Cylinder 1","Cylinder 2");
        xticks(0:60:720);
        xlim([0, 720]);
    end
    %}

    % New fourier coefficients from fourier approximations of cylinder torques
    Cn(:,[1,2,5])=0;
    Cn(:,2)=1/length(Mt)*fft(Mt(:,1));
    Cn(:,3)=1/length(Mt)*fft(Mt(:,2));

    for n = 1:24
        bn = [zeros(5,1) ; M\transpose(Cn(n+1,1:5))]; % Excitation Vector
        Fn = i*n*w*eye(10) - A; % Frequency Matrix
        gn = Fn\bn; % Frequency Response Vector
        bigTheta = 2*abs(gn);
        phi = angle(gn);
        for j=1:5
            % Store fourier factors for Eq. 22 summation in paper
            amp(rpm_num,n,j) = max(bigTheta(j)*cos((n/2)*alpha - phi(j)))*180/pi;
        end
    end
    theta(rpm_num,1:5) = sum(amp(rpm_num,:,:)); % Eq. 22, with j=1(1)4, max amplitudes for each element of crankshaft

    % Actuating dynamic torques
    T(rpm_num,1:2) = abs(theta(rpm_num,[2 5])-theta(rpm_num,[3 4])).*kt([2 4])*pi/180;
end

%
% plot amplitude of vibration vs rpm at end of crankshaft
figure;
hold on;
plot(rpm, theta(:,1), 'k-');
styles = {'r--' 'b--' 'g--' 'c--' 'm--'};
terms = [2 8 10 12 14]; % Fourier terms to graph
for n = 1:5
    plot(rpm, amp(:,terms(n),1), styles{1,n});
end
title("Torsional Vibration Amplitude at Crankshaft Sprockets");
ylabel("Amplitude (degrees)");
xlabel("Speed (rpm)");
xlim([800 6000]);
legend(["Total" string(terms/2)]);
hold off;
%}
%{
% plot amplitude of vibration vs rpm at flywheel (roughly equivalent without dampers)
figure;
hold on;
plot(rpm, theta(:,5), 'k-');
styles = {'r--' 'b--' 'g--' 'c--' 'm--'};
terms = [2 8 10 12 14]; % Fourier terms to graph
for n = 1:5
    plot(rpm, amp(:,terms(n),5), styles{1,n});
end
title("Torsional Vibration Amplitude at Flywheel");
ylabel("Amplitude (degrees)");
xlabel("Speed (rpm)");
xlim([800 6000]);
legend(["Total" string(terms/2)]);
hold off;
%}

%
% plot actuating torques vs rpm at ends of crankshaft
figure;
plot(rpm, T(:,1)/1000, 'b-', rpm, T(:,2)/1000, 'r-');
title("Actuating Torques on Crankshaft");
ylabel("Torque (kNm)");
xlabel("Speed (rpm)");
xlim([800 6000]);
legend("Sprocket End","Flywheel");
%}