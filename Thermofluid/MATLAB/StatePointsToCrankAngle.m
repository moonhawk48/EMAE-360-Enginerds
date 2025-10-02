
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % P1--------Intake Pressure
    % T1--------Intake Temp
    % cr--------Compression Ratio
    % Cp--------Heat Capacity Const Pressure
    % Cv--------Heat Capacity Const Vol
    % eta-------Combustion Efficiency
    % Q---------Gasoline Heating Value (J/kg)
    % OF--------Air to Fuel Ratio
    % D---------Displacement (m^3)
    % N---------Engine Speed (rev/min)
    % m_air-----Mass of air
libLoc = "C:\Users\Lucian\AppData\Roaming\CoolProp";
CA = (1:720)'; %crank angle array
C = 2; % Number of Cylinders
theta = deg2rad(CA);
c = (L/2)*(R+1-cos(theta)-(R^2-sin(theta).^2).^(1/2)); %clearance height state array
dV = pi*B^2*c/4;
V = Vc + dV;  

CA_spark = 361;
CA_b = 40;

T_inf = 300;
P_inf = 101300;
T = zeros(721,1);
P = zeros(721,1);
rho = zeros(721,1);
throttle = 1;
AFR = 10+(6)*throttle;
m_mix = 0;

%Fuel composition
x_oct = .76;
x_hept = .14;
x_eth = 0.1;
%air composition
x_N2 = 0.78;
x_O2 = 0.21;
x_Ar = 0.0093;
R = 287.052874; % Specific Gas Const

mm_air = getFluidProperty(libLoc,'M','T',T_inf,'P',P_inf,'Air');
mm_fuel = getFluidProperty(libLoc,'M','T',T_inf,'P',P_inf,'Octane')*x_oct + ...
    getFluidProperty(libLoc,'M','T',T_inf,'P',P_inf,'nHeptane')*x_hept+ ...
    getFluidProperty(libLoc,'M','T',T_inf,'P',P_inf,'Ethanol')*x_eth;
rho_mix = mm_mix*P_inf/(8.314*T_inf);

for i = 1:length(CA)
    disp(i);
    if i == 1
        %should be changed to combustion products
        T(i) = T_inf;
        P(i) = P_inf;
        rho(i) = getFluidProperty(libLoc,"D","T",T_inf,"P",P_inf,"Air");
    elseif i <= 180
         %needs valve characteristics
        rho_new = rho(2)*(1-(1/(AFR+1)))+rho_mix/(AFR+1);
        rho(i) = rho_new;
        T(i) = T(i-1);
        P(i) = P(i-1);
        dm_a = rho(i)*dV(i);
        m_mix = m_mix+dm_a;
        
    elseif i > 180 && i <= 360
        dCR = dV(i-1)/dV(180); %compression difference between crank angles
        rho(i) = rho(180)*dCR;
        P(i) = P(180)*(rho(i)/rho(i-1))^y;
        T(i) = P(i)/(rho(i)*R);

    elseif i > 360 && i <= 520
            dCR = dV(i)/dV(361) %compression difference between crank angles
            rho(i) = rho(359)*dCR

            if i >= CA_spark && i <= CA_b+CA_spark
                x_b = 1-exp(-5*((CA(i+1)-CA_spark)/CA_b)^3)
                reactants = x_b*m_mix*[x_N2,x_O2,x_Ar,x_oct,x_hept,x_eth];
                inputs = [1/AFR,rho(i),T(i),reactants]
                CEAout = CEAcall(inputs);
                T(i) = CEAout(1);
                P(i) = CEAout(2);
                rho(i) = CEAout(3);

            else
                P(i) = P(360)*(rho(i)/rho(i-1))^y;
                T(i) = P(i)/(rho(i)*R);

            end
    elseif i > 520 && i <= 720
            T = T(i-1);
            P = P(i-1);
            dm_a = rho(i)*dV(i);
            m_mix = m_mix-dm_a;

    end
    y = 1.38-0.2*exp(-900/T(i)); %Egnell function
end







    % Work
    % W = (Cv*(T3-T2) - Cv*(T4-T1))*(D*p1*(cr/(cr-1)));
    % %Wcyl = W/C;
    % 
    % % Power
    % P = 0.00134102*N*W/120;
    % 
    % % Specific Fuel Consumption
    % SFC = (C/W)*(m_mix/OF)*3.6*10^9; % In kg/J. If you want to convert to g/kW*Hr, multiply by 3600000000
    % 
    % 

