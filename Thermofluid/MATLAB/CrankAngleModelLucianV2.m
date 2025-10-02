
clear 
clc
libLoc = "C:\Users\Lucian\AppData\Roaming\CoolProp";

% Engine geometry
Vd_total = .00157;       % total displacement m^3
Ncyl = 2;                % number of cylinders 
Vd = Vd_total/Ncyl;      % displacement per cylinder
B = (4*Vd/pi)^(1/3);     % bore
S = B;                   % stroke (square design)
a = S/2;                 % crankshaft radiusdoub
L = (S+a)*1.1;           % conrod w/ 10% clearance
R = L/a;                 % Ratio of conrod to crank
CR = 9;                  % compression ratio
Vc = Vd/(CR-1);          % clearance volume/minimum volume m^3
A_piston = pi*(B/2)^2;    %piston head area

%valve geometry
Dv = 0.43*B;

Cd = 0.8; %discharge coefficient

% Operating conditions
N_init = 800;
throttle = 1; 
omega = 2*pi*N_init/60;  % in rad/s
mean_speed = 2*L*N_init;
CA = (1:720)'; %crank angle array
theta = deg2rad(CA);
c = (L/2)*(R+1-cos(theta)-(R^2-sin(theta).^2).^(1/2)); %clearance height state array
dV = pi*B^2*c/4;
V = Vc + dV; %Cylinder Volume State Array

%Fuel composition
x_oct = .76;
x_hept = .14;
x_eth = 0.1;
%air composition
x_N2 = 0.78;
x_O2 = 0.21;
x_Ar = 0.0093;

% Valve & injection timings in degrees.
%Cylinder pressure = atm around 50 degrees 
IVO = 45; IVC = 164;      
EVO = 540; EVC = 720;     
SOI = 0; EOI = 180;     
comb_start = 360;         
Comb_duration = 30;  



% Fluid properties
P_inf = 1e5;   
T_inf = 300; 

P_exh = 1.05e5; 
AFR = 14.7;       % stoichiometric air/fuel ratio              
R = 287;                 % J/kg-K for air
rho_inf = getFluidProperty(libLoc,'D','T',T_inf,'P',P_inf,'Air');
mm_air = getFluidProperty(libLoc,'M','T',T_inf,'P',P_inf,'Air');

mm_fuel = getFluidProperty(libLoc,'M','T',T_inf,'P',P_inf,'Octane')*x_oct + ...
    getFluidProperty(libLoc,'M','T',T_inf,'P',P_inf,'nHeptane')*x_hept+ ...
    getFluidProperty(libLoc,'M','T',T_inf,'P',P_inf,'Ethanol')*x_eth;

mm_mix = mm_fuel*1/(AFR+1)+mm_air*(1-(1/(AFR+1)));
rho_mix = mm_mix*P_inf/(8.314*T_inf);

cp_air = getFluidProperty(libLoc,'CPMASS','T',T_inf,'P',P_inf,'Air');
cv_air = getFluidProperty(libLoc,'CVMASS','T',T_inf,'P',P_inf,'Air');

y = cp_air/cv_air;
Qv = 44*10^6;


%Potentially use Wiebe function and Woschni correlation for heat transfer coefficient?

%model single cylinder (as a function

%model multiple cylinders working in tandem
%model multiple cylinders over time

%Wiebe function 9.3 connects crank angle to burn fraction to calculate how
%long the burn occurs over

%q_dot to the walls can be calculated by equation on page 676

%burn fraction can be calcualted by 9.23 but requires polytropic modelling

%equation 9.7 provides full energy balance for the combustion
%equation 9.8 provides more complicated but does not require polytropic
%constant

%time for crank degree
dt = 1/(2*pi*N_init/60);

for i = 1:length(CA)-1
    %kinematics
    th = deg2rad(i);
    x_piston = (a+L)-(a*cos(th)+sqrt(L^2-(a*sin(th))^2));
    w = 2*pi*N_init/60;
    vel_piston = w*(a*sin(th)+(sin(th)*cos(th)*a^2)/sqrt(L^2-(a*sin(th)))^2);
    %Egnell function
    y=1.38-0.2*exp(-900/T_inf);

    if i == 1
        %based on Initial state point analysis, update as needed
        P_cyl = 5e5; 
        T_cyl = 1500;
        m_air = 0;
    else
        %cylinder partial pressure and temperatre due to kinematics
            P_cyl = [P_cyl,P_cyl(i-1)*(V(i-1)/V(i))^y];
            T_cyl = [T_cyl,T_cyl(i-1)*(V(i-1)/V(i))^(y-1)];
    end


    %static exhaust gasses (which do not get expelled)
    %comp = CEA('problem','uv','equilibrium','f/a',AFR,'eq','rho,kg/m3',rho_mix,'reac', ...
    % 'ox','N2','wt%',x_N2*m_air,'t(k)',T_cyl,'ox', ...
    % 'O2','wt%',mO2,'t(k)',T_cyl, ...
    % 'ox','Ar','wt%',mAr,'t(k)',T_cyl, ...
    % 'fu','C8H18,isooctane','wt%',mOctane,'t(k)',T_cyl, ...
    % 'fu','C7H16,n-heptane','wt%',mHeptane,'t(k)',T_cyl, ...
    % 'fu','C2H5OH','wt%',mEthanol,'t(k)',T_cyl, ...
    % 'output','end');
    % 

    if i >= IVO && i <= IVC
        disp(i)
            if i <= IVO + 0.1*(IVC-IVO) || i >= IVO + 0.9*(IVC-IVO)
                state = 1;
            else
                state = 2;
            end
        if P_cyl(i) <= P_inf
            [mdot_a,choked] = Valve(state,B,T_inf,P_inf,P_cyl(i),"Air",libLoc)
            dm_a = mdot_a*dt*throttle;
            dm_f = dm_a/AFR;
            dm_mix = dm_f+dm_a;
            m_air = m_air + dm_a;
            %cylinder partial pressure + air partial pressure
            P_cyl(i) = P_cyl(i)+(dm_mix/mm_mix)*8.314*T_inf/dV(i);
        else %backflow through intake (very bad)
            [mdot_a,choked] = Valve(state,B,T_inf,P_cyl(i),P_inf,"Air",libLoc)
            dm_a = mdot_a*dt*throttle;
            dm_f = dm_a/AFR;
            dm_mix = dm_f+dm_a;
            m_air = m_air - dm_a
            %cylinder partial pressure + air partial pressure
            P_cyl(i) = P_cyl(i-1)-(dm_mix/mm_mix)*8.314*T_inf/dV(i);
            if P_cyl(i) <= 0
               P_cyl(i) = 0;
               m_air = m_air + dm_a
            end
        end
    end
end
plot(P_cyl)

%piston moves down, creates pressure differential which drives air via
%bernoulli's + flow restrictions
%throttle limits mass flow rate of air
%injectors mix fuel w air
%intake valve 