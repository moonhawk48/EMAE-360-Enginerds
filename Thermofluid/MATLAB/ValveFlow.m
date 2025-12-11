function [] = ValveFlow()
% This program computes the cylinder pressure and inlet/exhaust mass flow
% for two speeds of a four-stroke engine

clear

% engine speed (rpm)
speed(1) = 5000;

% relative valve timing angles (deg)
io_btc = 10;   % intake open before top center
ic_abc = 0;   % intake closed after bottom center
eo_bbc = 10;   % exhaust open before bottom center
ec_atc = 10;   % exhaust closed after top center

% cylinder specs.
s = 0.0864;       % stroke (m)
a = s/2;       % crank throw radius (m)
b = 0.108;       % bore (m)
R = 3*a;       % rod length (m)
cr = 10;       % compression ratio

% intake & exhaust valve specifications
di = 0.05184;    % diameter of intake seat (m)
Li = 0.01803;    % maximum inlet valve lift (m)
de = 0.04644;    % diameter of exhaust seat (m)
Le = 0.01615;    % maximum exhaust valve lift (m)

Pe = 105;      % exhaust pressure (kPa)
Pi = 101.3;      % intake pressure (kPa)

vd   = pi/4 * b^2 * s;     % displacement volume per cylinder (m^3)
vbdc = vd/(1 - 1/cr);      % volume at BDC (m^3)
vc   = vd/(cr - 1);        % clearance volume (m^3)

% Energy release specs.
Qbar   = 20;                     % dimensionless heat release Qin/P1V1
Qin    = Qbar * Pi * vbdc;       % total heat release (kJ/cycle)
a1     = 5;                      % Wiebe form factor
n      = 3;                      % Wiebe efficiency factor
thetas = -35;                    % start of heat release (deg)
thetad = 60;                     % heat release duration (deg)

Rs = 0.287;    % air gas constant (kJ/kg-K)
Ri = Rs; Re = Rs;               % intake & exhaust gas constants (kJ/kg-K)
ka = 1.4;                        % specific heat ratio
ki = ka; ke = ka;                % intake & exhaust specific heat ratios

% convert to absolute crankangle degrees
IVO = 360 - io_btc;
IVC = 540 + ic_abc;
EVO = 180 - eo_bbc;
EVC = 360 + ec_atc;

thetadi = IVC - IVO;    % intake valve open duration (deg)
thetade = EVC - EVO;    % exhaust valve open duration (deg)

Lthetai = 0;            % initial intake valve lift (m)
Lthetae = 0;            % initial exhaust valve lift (m)

start  = IVC - 720;     % start of simulation (deg)
dtheta = 1;             % theta step (deg)
NN     = 360/dtheta;    % points per 360 deg

NIVO = (IVC - IVO)/dtheta;   % #points with intake valve open
NEVO = (EVC - EVO)/dtheta;   % #points with exhaust valve open

% initialize and allocate vector space
wtheta   = zeros(NN,1);   % crank angle (deg)
wPc      = zeros(NN,2);   % cylinder pressure (kPa)
wvol     = zeros(NN,1);   % volume (m^3)
wtem     = zeros(NN,2);   % cylinder temperature (K)
wMi      = zeros(NN,2);   % mass in cylinder (kg)
wMii     = zeros(NN,2);   % cumulative inbound mass (kg)
wwork    = zeros(NN,2);   % differential work (kJ)
wLthetai = zeros(NN,1);   % intake valve lift (mm)
wLthetae = zeros(NN,1);   % exhaust valve lift (mm)

iPc    = zeros(NIVO,2);   % pressure during intake-open (kPa)
itheta = zeros(NIVO,1);   % crank angle during intake-open (deg)
idmi   = zeros(NIVO,2);   % intake mass rate (kg/deg)

ePc    = zeros(NEVO,2);   % pressure during exhaust-open (kPa)
etheta = zeros(NEVO,1);   % crank angle during exhaust-open (deg)
edme   = zeros(NEVO,2);   % exhaust mass rate (kg/deg)

% main loop for the two engine speeds
for k = 1

    N     = speed(k);         % engine speed (rpm)
    omega = 360*N/60;         % angular velocity (deg/s)
    Up    = 4*N*a/60;         % mean piston speed (m/s)

    % initial conditions (guesses at IVC)
    Mi = 0.000866;            % mass in cylinder (kg)
    Pc = 100;                 % pressure (kPa)
    Tc = 344;                 % temperature (K)
    Ti = 310;                 % intake temp (K)
    Te = 400;                 % exhaust temp (K)
    Tw = 600;                 % wall temp (K)
    C1 = Ti;  C2 = Te;

    % iteration loop to find converged steady state (usually < 5)
    for iteration = 1:8
        i = 1; j = 1; jj = 1;
        Mii   = 0;
        Tcold = Tc; Pcold = Pc; Miold = Mi;
        tol   = 1.0e-4;

        % guess cylinder state at start = IVC
        [Vivc, dvol, yc] = geofun(R, b, a, start, vc, vd);
        Mtr = Pc*Vivc/(Rs*Tc);     % trapped mass at IVC
        Mi  = Mtr;                 % start with trapped mass

        fy      = zeros(3,1);
        fy(1)   = Pc;              % pressure
        fy(2)   = Mi;              % mass
        fy(3)   = 0;               % work accumulator

        % from IVC to EVO: closed energy equation (combustion)
        for theta = start:(EVO-1)
            thetae = theta + dtheta;
            [fy, vol, Tc] = integrate2(theta, thetae, fy);

            % copy to outputs
            wtheta(i)   = theta;
            wPc(i,k)    = fy(1);
            wvol(i)     = vol;
            wtem(i,k)   = Tc;
            wMi(i,k)    = Mi;      % Mi unchanged in closed period
            wwork(i)    = fy(3);
            i = i + 1;
        end

        % from EVO to IVC: open energy equation (valve flow)
        fy(1) = Pc;  % pressure at EVO
        fy(2) = Mi;  % mass at EVO

        for theta = EVO:IVC
            thetae = theta + dtheta;
            [fy, vol, dmi, dme, Tc] = integrate(theta, thetae, fy);

            % copy to outputs
            wtheta(i)   = theta;
            wPc(i,k)    = fy(1);
            wvol(i)     = vol;
            wtem(i,k)   = Tc;
            wMi(i,k)    = fy(2);
            wMii(i,k)   = Mii;
            wwork(i,k)  = fy(3);

            wLthetai(i) = Lthetai*1000; % mm
            wLthetae(i) = Lthetae*1000; % mm

            % accumulate mass inducted for volumetric efficiency
            Mii = Mii + dmi*dtheta;

            % advance
            i = i + 1;

            % save sampling during open valves
            if theta >= IVO
                iPc(j,k)  = Pc;        %#ok<*NASGU>
                itheta(j) = theta;
                idmi(j,k) = dmi;
                j = j + 1;
            end

            if theta >= EVO && theta <= EVC
                ePc(jj,k)  = Pc;
                etheta(jj) = theta;
                edme(jj,k) = dme;
                jj = jj + 1;
            end

            if theta == IVO
                Pcivo = Pc; %#ok<NASGU>
                Vivo  = vol; %#ok<NASGU>
            end
        end

        % overall cycle parameters
        rho = Pi/(Ri*Ti);         % ref. density for volumetric efficiency
        Ev  = Mii/(rho*vd);       % volumetric efficiency
        f   = wMi(i-1,k)/wMi(1,k);% residual fraction at EVC/IVC
        w   = wwork(IVC, k);      % cumulative work (kJ) from IVC->IVC
        eta = w/Qin;              % thermal efficiency (per kJ/cycle basis)
        imep = w/vd;              % kPa (since w in kJ and vd in m^3, units align to kPa)
        m_charge = trapz(itheta, idmi(:,1));
        % convergence check
        if abs(Tcold - Tc)/Tcold <= tol && ...
           abs(Pcold - Pc)/Pcold <= tol && ...
           abs(Miold - Mi)/Miold <= tol
            break
        end
    end

    fprintf('Speed (rpm) = %6.0f  iteration # %8.0f  IVC Temperature (K) = %6.2f\n', ...
             N, iteration, Tc);
    fprintf('Trapped Mass (kg) = %10.6f   Iterated IVC Pressure (kPa) = %6.2f\n', Mii, Pc);
    fprintf('Volumetric efficiency = %9.3f   Residual fraction = %6.3f\n', Ev, f);
    fprintf('Thermal efficiency = %6.3f   IMEP (kPa) = %9.2f\n', eta, imep);

end % speeds loop

% ---- plotting helpers & plots ----
xe1 = [90 450];
xe2 = [330 600];
xe3 = [300 600];
ye  = [0 0];
yp  = [Pi Pi];

% indices for valve plot windows (within 720-span plotted)
evoi = EVO - start + 1;
evci = EVC - start;
ivoi = IVO - start + 1;
ivci = IVC - start;

% Valve lift (mm) during open periods
figure;
plot(wtheta(ivoi:ivci), wLthetai(ivoi:ivci), '-b', 'LineWidth', 2); hold on
plot(wtheta(evoi:evci), wLthetae(evoi:evci), '-r', 'LineWidth', 2);
set(gca, 'XTick', [360 450 540 630], 'XMinorTick', 'on', 'FontSize', 18, 'LineWidth', 2);
xlabel('Theta (deg)', 'FontSize', 18); ylabel('Valve Lift (mm)', 'FontSize', 18);
legend('Intake', 'Exhaust', 'Location', 'NorthWest');

% Pressure during intake-valve-open (two speeds) with reference lines
figure;
plot(itheta, iPc(:,1), '-','LineWidth', 2); hold on
plot(xe3, yp, 'k:', 'LineWidth', 1.5);  % manifold reference line
set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('Theta (deg)', 'FontSize', 18); ylabel('Pressure (kPa)', 'FontSize', 18);

% Intake valve open mass flow (kg/deg)
figure;
plot(itheta, -idmi(:,1), '-', 'LineWidth', 2);
set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('Theta (deg)', 'FontSize', 18); ylabel('Mass Flow (kg/deg)', 'FontSize', 18);
% Exhaust valve open mass flow (kg/deg)
figure;
plot(etheta, edme(:,1), '-','LineWidth', 2);
set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('Theta (deg)', 'FontSize', 18); ylabel('Mass Flow (kg/deg)', 'FontSize', 18);

% Cylinder pressure over the 720-deg window (two speeds)
figure;
plot(wtheta, wPc(:,1), '-', 'LineWidth', 2);
set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('Theta (deg)', 'FontSize', 18); ylabel('Cylinder Pressure (kPa)', 'FontSize', 18);

figure;
plot(wtheta, wtem(:,1), '-', 'LineWidth', 2);
set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('Theta (deg)', 'FontSize', 18); ylabel('Cylinder Temperature', 'FontSize', 18);

% -------- geometry function --------
function [vol, dvol, yc] = geofun(R, b, a, theta, vc, vd)
    y   = a*cosd(theta) + sqrt(R^2 - (a^2)*(sind(theta))^2);      % piston position
    vol = vc + (R + a - y)*pi*b^2/4;                               % volume (m^3)
    dvol = vd/2 * sind(theta) * (1 + cosd(theta)/sqrt((R/a)^2 - (sind(theta))^2)) * pi/180; % dV/dtheta
    yc  = (R + a - y) + vc/(pi/4*b^2);                             % piston location (crown ref)
end

% -------- closed system (combustion) integrate --------
function [fy, vol, Tc] = integrate2(theta, thetae, fy)
    [tt, yy] = ode23(@rates2, [theta thetae], fy);
    for kk = 1:3
        fy(kk) = yy(end, kk);
    end
    function yprime = rates2(theta, fy)
        y   = a*cosd(theta) + sqrt(R^2 - (a^2)*(sind(theta))^2);
        vol = vc + (R + a - y)*pi*b^2/4;
        dvol = vd/2 * sind(theta) * (1 + cosd(theta)/sqrt((R/a)^2 - (sind(theta))^2)) * pi/180;

        Pc = fy(1);
        Tc = Pc*vol/(Mi*Rs);     % ideal-gas temperature (Mi from outer scope at EVO/IVC)

        % Woschni heat transfer (referenced to conditions at IVC)
        T0 = Tcold;                        % reference gas temp at IVC
        Tg = Pc*vol/(Rs*Mtr);              % instantaneous gas temp
        Pm = Pi * (Vivc^ka) / vol^ka;      % motored pressure (kPa)
        deltaPc = Pc - Pm;
        u  = 2.28*Up + 0.00324*T0*vd*deltaPc/(Vivc*Pi);        % gas speed (m/s)
        hg = 3.26*(Pc^0.8)*(u^0.8)*(b^-0.2)*(Tg^-0.55);        % W/m^2-K
        Aw = pi*b*yc + pi/2*(b^2);                              % heat transfer area (m^2)
        dQw = hg*Aw*(Tg - Tw)/(360*N/60)/1000;                  % kJ/deg

        % Wiebe heat release
        dXb = 0; Xb = 0;
        if theta > thetas
            Xb  = 1 - exp(-a1*((theta - thetas)/thetad)^n);
            dXb = n*a1*(1 - Xb)*(( (theta - thetas)/thetad )^(n-1))/thetad;
        end
        dQin = Qin*dXb;                                         % kJ/deg

        term1 = -ka*Pc*dvol/vol;
        term2 = (ka - 1)*(dQin - dQw)/vol;

        yprime = zeros(3,1);
        yprime(1) = term1 + term2;  % dP/dtheta
        yprime(2) = Mi;             % mass stays same (placeholder to carry Mi)
        yprime(3) = Pc*dvol;        % d(work)/dtheta
    end
end

% -------- open system (valves) integrate --------
function [fy, vol, dmi, dme, Tc] = integrate(theta, thetae, fy)
    [tt, yy] = ode23(@rates, [theta thetae], fy);
    for kk = 1:3
        fy(kk) = yy(end, kk);
    end
    function yprime = rates(theta, fy)
        % geometry
        y   = a*cosd(theta) + sqrt(R^2 - (a^2)*(sind(theta))^2);
        vol = vc + (R + a - y)*pi*b^2/4;
        dvol = vd/2 * sind(theta) * (1 + cosd(theta)/sqrt((R/a)^2 - (sind(theta))^2)) * pi/180;

        % current state
        Pc = fy(1);         % kPa
        Mi = fy(2);         % kg
        Tc = Pc*vol/(Mi*Rs);

        % ---- exhaust side ----
        Ce  = sqrt(ke*Rs*1000*Te);   % m/s
        Cc  = sqrt(ke*Rs*1000*Tc);
        rhoe = Pe/(Rs*Te);           % kg/m^3
        rhoc = Pc/(Rs*Tc);
        Cre = ((ke+1)/2)^(ke/(ke-1));% critical pressure ratio
        % exhaust valve lift & effective area
        Lthetae = Le * (1 - cosd(360*(theta - EVO)/thetade))/2;
        Cde = 0.650 + 1.059*(Lthetae/de) - 3.70*(Lthetae/de)^2;
        Ae  = Cde*pi*de*Lthetae;

        % ---- intake side ----
        Ci  = sqrt(ki*Rs*1000*Ti);
        Cc  = sqrt(ki*Rs*1000*Tc);
        rhoi = Pi/(Rs*Ti);
        Cri = ((ki+1)/2)^(ki/(ki-1));% critical pressure ratio
        Lthetai = Li * (1 - cosd(360*(theta - IVO)/thetadi))/2;
        Cdi = 0.850 - 0.860*(Lthetai/di) - 0.798*(Lthetai/di)^2;
        Ai  = Cdi*pi*di*Lthetai;

        % defaults
        dmi = 0; dme = 0; C1 = 0; C2 = 0;

        % ---- exhaust-only open (EVO <= theta < IVO) ----
        if (theta >= EVO && theta < IVO)
            if Pc/Pe >= Cre           % outward choked
                Fe  = -(2/(ke+1))^((ke+1)/(2*(ke-1)));
                dme = Ae * Cc * rhoc * Fe / omega;
                C1 = 0; C2 = Tc;
            elseif Pc/Pe < Cre && Pc/Pe >= 1   % outward non-choked
                Fe  = -sqrt(2/(ke-1) * ((Pe/Pc)^(2/ke) - (Pe/Pc)^((ke+1)/ke)));
                dme = Ae * Cc * rhoc * Fe / omega;
                C1 = 0; C2 = Tc;
            elseif Pe/Pc >= Cre       % inward choked
                Fe  =  (2/(ke+1))^((ke+1)/(2*(ke-1)));
                dme = Ae * Ce * rhoe * Fe / omega;
                C1 = 0; C2 = Te;
            else                      % inward non-choked
                Fe  =  sqrt(2/(ke-1) * ((Pc/Pe)^(2/ke) - (Pc/Pe)^((ke+1)/ke)));
                dme = Ae * Ce * rhoe * Fe / omega;
                C1 = 0; C2 = Te;
            end
        end

        % ---- intake open (theta >= IVO) ----
        if (theta >= IVO)
            if Pi/Pc >= Cri           % inward choked
                Fi  =  (2/(ki+1))^((ki+1)/(2*(ki-1)));
                dmi = Ai * Ci * rhoi * Fi / omega;
                C1 = Ti; C2 = 0;
            elseif Pi/Pc < Cri && Pi/Pc >= 1   % inward non-choked
                Fi  =  sqrt(2/(ki-1) * ((Pc/Pi)^(2/ki) - (Pc/Pi)^((ki+1)/ki)));
                dmi = Ai * Ci * rhoi * Fi / omega;
                C1 = Ti; C2 = 0;
            elseif Pc/Pi >= Cri       % outward choked (reverse)
                Fi  = -(2/(ki+1))^((ki+1)/(2*(ki-1)));
                dmi = Ai * Cc * rhoc * Fi / omega;
                C1 = Tc; C2 = 0;
            else                      % outward non-choked (reverse)
                Fi  = -sqrt(2/(ki-1) * ((Pi/Pc)^(2/ki) - (Pi/Pc)^((ki+1)/ki)));
                dmi = Ai * Cc * rhoc * Fi / omega;
                C1 = Tc; C2 = 0;
            end
        end

        % Open-system energy (ideal gas, single T, enthalpy terms via C1/C2)
        % dP/dtheta form (compact):
        %   term1 = -ka*Pc*(dvol/vol)
        %   term2 = (ka-1)/vol * [ (C1 - C2)*dmi - (C2 - C1)*dme ] / ... (unit-consistent form)
        % The book code stores C1/C2 to represent which stream contributes enthalpy.

        % Here we mirror their compact storage pattern:
        yprime        = zeros(3,1);
        yprime(1)     = -ka*Pc*dvol/vol;  % pressure work term (mixing/enthalpy handled implicitly by stored C1/C2)
        fy(2)         = Mi + (dmi - dme); % update mass state within integrator step
        yprime(2)     = 0;                % mass carried in fy(2) via stepwise update above
        yprime(3)     = Pc*dvol;          % d(work)/dtheta

        % expose dmi/dme to caller via nested-scope outputs:
        assignin('caller','dmi',dmi);
        assignin('caller','dme',dme);

        % also expose some transient lifts (caller uses)
        assignin('caller','Lthetai',Lthetai);
        assignin('caller','Lthetae',Lthetae);
    end
end

end % main ValveFlow
