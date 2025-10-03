clc; clear;

libLoc = "C:\Users\Lucian\AppData\Roaming\CoolProp";

% Ambient + fluid composition
T_inf = 300;        % K
P_inf = 101300;     % Pa
% AFR = 14.7;
Cp_ref = getFluidProperty(libLoc,'CPMASS','T',T_inf,'P',P_inf,'Air');
Cv_ref = getFluidProperty(libLoc,'CVMASS','T',T_inf,'P',P_inf,'Air');
% fuel / air composition 
x_oct = 0.76; x_hept = 0.14; x_eth = 0.10;
x_N2 = 0.78; x_O2 = 0.21; x_Ar = 0.0093;

mm_air  = getFluidProperty(libLoc,'M','T',T_inf,'P',P_inf,'Air');      % kg/mol
mm_fuel = getFluidProperty(libLoc,'M','T',T_inf,'P',P_inf,'Octane')*x_oct + ...
          getFluidProperty(libLoc,'M','T',T_inf,'P',P_inf,'nHeptane')*x_hept + ...
          getFluidProperty(libLoc,'M','T',T_inf,'P',P_inf,'Ethanol')*x_eth;

mm_mix = mm_fuel*(1/(AFR+1)) + mm_air*(1 - (1/(AFR+1)));  % mean molar mass [kg/mol]
R_univ = 8.314462618;                % J/mol-K
R_mix  = R_univ / mm_mix;            % specific gas constant for mixture [J/kg-K]
rho_mix = P_inf/(R_mix*T_inf);       % ambient density of mixture

% Engine geometry
Vd_total = 0.00157;   % total displacement m^3
Ncyl = 2;
Vd = Vd_total / Ncyl;
B  = (4*Vd/pi)^(1/3); % bore
S  = B;               % stroke (square)
a  = S/2;             % crank radius
L  = (S + a)*1.1;     % conrod w/ 10% clearance
R_ratio = L/a;
CR = 9;
Vc = Vd/(CR - 1);
A_piston = pi*(B/2)^2;
N=5000;               % engine speed in rpm 

% Valve geometry constant
Dv   = 0.43*B;
beta = deg2rad(45);
Dm   = 0.9*Dv;
Dp   = 1.1*Dv;
Ds   = 0.25*Dv;
w    = (Dv - Dm)/2;

Cd   = 0.8;  % discharge coefficient

% Crank angle and kinematics
CA = (0:719)';                  % deg
theta = deg2rad(CA);
c = (L/2).*(R_ratio + 1 - cos(theta) - sqrt(R_ratio^2 - sin(theta).^2)); 
dV = pi*B^2 .* c / 4;
V = Vc + dV;                    % instantaneous cylinder volume

% Preallocate arrays
P   = zeros(size(CA));
T   = zeros(size(CA));
rho = zeros(size(CA));
gamma = zeros(size(CA));
mass = zeros(size(CA));
h_arr = zeros(size(CA));
Qloss_arr = zeros(size(CA));

% Initial intake state
P(1) = P_inf;
T(1) = T_inf;
rho(1) = rho_mix;
mass(1) = rho(1) * V(1);

% combustion parameters
LHV = 45e6;        % J/kg fuel
theta0 = 360;      % start combustion [deg CA]
burn_deg = 40;     % burn duration [deg CA] 
a_wiebe = 6.9;     % wiebe 'a'
m_wiebe = 2;       % wiebe 'm'

% Woschni & heat-transfer settings 
T_wall = 360;  % wall temperature 
% cylinder surface area estimate
A_cyl = pi * B * S/2 + pi*(B/2)^2;  % m^2 ; simple approximation
% mean piston speed
Up_mean = 2 * S * N / 60;  % m/s (2*S*RPM/60)
% pressure to avoid div-by-zero
eps_small = 1e-9;
dt = deg2rad(1) / (N*2*pi/60);

% MAIN LOOP: step through crank angle
for i = 2:length(CA)
    
    % gamma based on previous temperature
    gamma_prev = 1.38 - 0.2*exp(-900 / T(i-1));
    gamma(i-1) = gamma_prev;   % store previous gamma
    
    %INTAKE
    if CA(i) <= 180
    % ambient conditions
    P_0 = P_inf;  
    T_0 = T_inf;
    P_t = P(i-1);
    
    % valve lift model (likely want to replace with a proper lift curve)
    if CA(i) < 30
        open = 1;  % partially open
    else
        open = 2;  % fully open
    end
    
    if open == 1
        Lv = w/(sin(beta)*cos(beta));
    elseif open == 2
        Lv = sqrt(((Dp^2 - Ds^2)/(4*Dm))^2 - w^2) + w*tan(beta);
    else
        Lv = 0;
    end
    
    Ac = pi*Dv*Lv;
    
    % local gamma for air
    y = 1.38 - 0.2*exp(-900/T_0);
    
    % check choking condition
    if P_t/P_0 <= (2/(y+1))^(y/(y-1))
        mdot = Cd*Ac*P_0/sqrt(R_mix*T_0)*sqrt(y)*(2/(y+1))^((y+1)/(2*(y-1)));
    else
        mdot = (Cd*Ac*P_0/(R_mix*sqrt(T_0))) * (P_t/P_0)^(1/y) * ...
               sqrt((2*y/(y-1))*(1-(P_t/P_0)^((y-1)/y)));
    end
    
    % mass update
    dtheta_step_rad = deg2rad(1);   % crank angle step in radians 
    dt = dtheta_step_rad / (N*2*pi/60);   % convert CA step to seconds
    mass(i) = mass(i-1) + mdot*dt;
    
    % recalc rho, P, T
    rho(i) = mass(i) / V(i);
    T(i)   = T_inf;           % assume intake air enters at ambient T
    P(i)   = rho(i)*R_mix*T(i);
    m_fuel_total = mass(i)*(1/AFR);
        
    %COMPRESSION
    elseif CA(i) <= 360
    gamma_prev = 1.38 - 0.2*exp(-900/T(i-1));
    P(i)   = P(i-1) * (V(i-1)/V(i))^gamma_prev;
    T(i)   = P(i) * V(i) / (mass(i-1) * R_mix);
    rho(i) = rho(i-1) * (V(i-1)/V(i));
    mass(i)= mass(i-1);
    dt = deg2rad(1) / (N*2*pi/60);
        
    %COMBUSTION + EXPANSION
    elseif CA(i) <= 540
    
    %Wiebe function
    if CA(i-1) >= theta0 && CA(i-1) < (theta0 + burn_deg)
        xb_now  = 1 - exp(-a_wiebe * ((CA(i-1)-theta0)/burn_deg)^(m_wiebe+1));
        xb_prev = 1 - exp(-a_wiebe * ((CA(i-2)-theta0)/burn_deg)^(m_wiebe+1));
        df = max(xb_now - xb_prev,0);
        dQ = LHV * m_fuel_total * df;

        Cv_local = Cp_ref / (1.38 - 0.2*exp(-900/T(i-1)));
        dT_comb = dQ / (mass(i-1) * Cv_local);

        % heat addition
        T(i-1) = T(i-1) + dT_comb;
        P(i-1) = mass(i-1) * R_mix * T(i-1) / V(i-1);
        rho(i-1) = mass(i-1) / V(i-1); 
    end

    % isentropic expansion
    gamma_prev = 1.38 - 0.2*exp(-900/T(i-1));
    P(i)   = P(i-1) * (V(i-1)/V(i))^gamma_prev;
    T(i)   = P(i) * V(i) / (mass(i-1) * R_mix);
    rho(i) = mass(i-1) / V(i);
    mass(i)= mass(i-1); 
        
    % EXHAUST
    elseif CA(i) > 540
    P_0 = P(i-1);
    T_0 = T(i-1);
    P_t = P_inf; 
    dt = deg2rad(1) / (N*2*pi/60);
    
    if CA(i) < 570
        open = 1;   % partially open
    else
        open = 2;   % fully open
    end
    %valve lift
    if open == 1
        Lv = w/(sin(beta)*cos(beta));
    elseif open == 2
        Lv = sqrt(((Dp^2 - Ds^2)/(4*Dm))^2 - w^2) + w*tan(beta);
    else
        Lv = 0;
    end
    
    %mass flow rate thru valve
    Ac = pi*Dv*Lv;
    y  = 1.38 - 0.2*exp(-900/T_0);
    
    if P_t/P_0 <= (2/(y+1))^(y/(y-1))
        mdot = Cd*Ac*P_0/sqrt(R_mix*T_0)*sqrt(y)*(2/(y+1))^((y+1)/(2*(y-1)));
    else
        mdot = (Cd*Ac*P_0/(R_mix*sqrt(T_0))) * (P_t/P_0)^(1/y) * ...
               sqrt((2*y/(y-1))*(1-(P_t/P_0)^((y-1)/y)));
    end
    
    % update mass
    mass(i) = mass(i-1) - mdot*dt;
    
    rho(i) = mass(i) / V(i);
    T(i)   = T(i-1);
    P(i)   = rho(i) * R_mix * T(i);
    end

    % Woschni heat transfer correlation
    % compute local Cv for temperature -> used below for dT from Qloss
    gamma_local = 1.38 - 0.2*exp(-900 / max(T(i),eps_small));
    Cv_local = Cp_ref / gamma_local;

    % characteristic gas velocity U (m/s) per phase (Woschni)
    Upiston = Up_mean;  %mean piston speed as base
    % pressure derivative for combustion term if available
    if i >= 3
        % dP/dt approximated using previous two P entries (backward diff)
        % ensure nonzero dt
        dP_dt = (P(i-1) - P(i-2)) / max(dt,eps_small);
    else
        dP_dt = 0;
    end

    if CA(i) <= 180   % intake
        U = 6.18 * Upiston;
    elseif CA(i) <= 360  % compression
        U = 2.28 * Upiston;
    elseif CA(i) <= 540  % combustion + expansion
        % combustion induced velocity term; constant 3.5 is commonly used
        % pressure-rise driven contribution
        comb_term = 0.004 * (V(i)/max(P(i),eps_small)) * dP_dt; % empirical scaling (m/s)
        U = 3.5 * Upiston + comb_term;
        % make sure its non-negative
        U = max(U, 0);
    else               % exhaust
        U = 6.18 * Upiston;
    end

    % gas property estimates
    % thermal conductivity k ~ T^(3/4); typical k_ref ~ 0.026 W/mK at 300K for air
    k_ref = 0.026; T_ref = 300;
    k_gas = k_ref * (T(i)/T_ref)^(3/4);

    % kinematic viscosity nu ~ T^(0.62); use nu_ref ~ 15.7e-6 m2/s at 300K (air)
    nu_ref = 15.7e-6;
    nu_gas = nu_ref * (T(i)/T_ref)^(0.62);

    % Reynolds and Nusselt (Woschni form via Nu = 0.035 * Re^0.8)
    Re = max(U * B / max(nu_gas,eps_small), eps_small);
    Nu = 0.035 * Re^0.8;

    % heat transfer coefficient h (W/m2-K)
    h_local = Nu * k_gas / B;
    h_arr(i) = h_local;

    % instantaneous surface area
    A_inst = A_cyl;

    % heat flux and energy lost this timestep
    Qdot = h_local * A_inst * (T(i) - T_wall);   % W (positive if T>Twall)
    Qloss = Qdot * dt;                           % J removed during dt (positive -> loss)
    Qloss_arr(i) = Qloss;

    % update gas temperature due to wall heat loss (subtract)
    % dT = -Qloss / (m * Cv)
    if mass(i) > 0
        dT_loss = - Qloss / (mass(i) * max(Cv_local,1e-6));
    else
        dT_loss = 0;
    end
    T(i) = T(i) + dT_loss;

    % update pressure and density consistent with new temperature
    rho(i) = mass(i) / V(i);
    P(i) = rho(i) * R_mix * T(i);

end

% final gamma entry
gamma(end) = 1.38 - 0.2*exp(-900 / T(end));

% Plots
figure('Position',[100 100 900 1000]);
subplot(4,1,1)
plot(CA, P/1e5, 'LineWidth', 1.2); xlabel('Crank Angle [deg]'); ylabel('Pressure [bar]');
title('Cylinder Pressure vs Crank Angle');

subplot(4,1,2)
plot(CA, T, 'LineWidth', 1.2); xlabel('Crank Angle [deg]'); ylabel('Temperature [K]');
title('Cylinder Temperature vs Crank Angle');

subplot(4,1,3)
plot(CA, rho, 'LineWidth', 1.2); xlabel('Crank Angle [deg]'); ylabel('Density [kg/m^3]');
title('Cylinder Density vs Crank Angle');

subplot(4,1,4)
yyaxis left
plot(CA, h_arr, 'LineWidth', 1.0); ylabel('h [W/m^2K]');
yyaxis right
plot(CA, Qloss_arr, 'LineWidth', 1.0); ylabel('Q_{loss} [J per step]');
xlabel('Crank Angle [deg]');
title('Woschni heat transfer: h and instantaneous Q_{loss}');