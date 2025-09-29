function [OttoCycleData] = Otto_V4_In_Sheet_State_Points_w_OF_and_CR(Tin, Pin, cr, Cp, Cv, eta, Q, OF, D, N)
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
    
    T1 = Tin;
    P1 = Pin;
    R = 287.052874; % Specific Gas Const
    g = Cp/Cv; % Heat Capacity Ratio
    C = 2; % Number of Cylinders

    % Stage 1: Intake
    p1 = P1/(R*T1); % Density

    % Stage 2: Compression
    p2 = p1*cr;
    P2 = P1*(p2/p1)^g;
    T2 = T1*(P2/P1)^((g-1)/g);

    % Stage 3: Combustion
    T3 = T2 + eta*(Q/(Cp*OF));
    p3 = p2;
    P3 = p3*R*T3;

    % Stage 4: Expansion
    p4 = p1;
    P4 = P3*(p4/p3)^g;
    T4 = P4/(p4*R);

    % Work
    W = (Cv*(T3-T2) - Cv*(T4-T1))*(D*p1*(cr/(cr-1)));
    %Wcyl = W/C;

    % Power
    P = 0.00134102*N*W/120;

    % Specific Fuel Consumption
    SFC = (3.6*10^9)/(OF*(Cv*(T3-T2) - Cv*(T4-T1))); % In kg/J. If you want to convert to g/kW*Hr, multiply by 3600000000

    %{
    DataTable = [   "Stage 1, Intake: "         P1 T1   p1
                    "Stage 2, Compression: "    P2 T2   p2
                    "Stage 3, Combustion: "     P3 T3   p3
                    "Stage 4, Expansion: "      P4 T4   p4
                    "Work: "                    W  Wcyl "N/A"
                    "Power, SFC: "              P  SFC  "N/A"]; % Aesthetically pleasing if debugging is necessary
    %}

    OttoCycleData = [P1; T1; p1; P2; T2; p2; P3; T3; p3; P4; T4; p4; W; P; SFC];

end