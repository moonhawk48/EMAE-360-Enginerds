function [Cp_mix, Cv_mix, h_mix] = mixProps(T, P, AFR, varargin)
%Inputs: mixture temp, pressure, AFR
%Optional inputs: mass of fuel air mixture, residual mass

FuelID      = 2;        % gasoline
cpLib       = '';       % CoolProp library path (if your wrapper needs it)
FS_gasoline = 0.06548;  % stoichiometric fuel/air mass ratio (fuel/air)
AFRst       = 1/FS_gasoline;

 P_kPa = P/1000;
 P_Pa  = P;


phi = AFRst / AFR;
[~, h_res, ~, ~, ~, R_res, Cp_res, ~, ~, ~, ~] = farg(T, P_kPa, phi, 1, FuelID);

Cp_air_J = getFluidProperty(cpLib, 'CPMASS','T',T,'P',P_Pa);
h_air_J  = getFluidProperty(cpLib, 'HMASS','T',T,'P',P_Pa);
Cp_air = Cp_air_J/1000;  % kJ/kg-K
h_air  = h_air_J /1000;  % kJ/kg
R_air  = 0.287042;       % kJ/kg-K (ideal-gas air)

% fuel() gives cp_fuel and h_fuel in *non-dimensional* (cp/R, h/RT) form in your tables.

[~, ~, ~, ~, h_fuel_nd, ~, cp_fuel_nd, MW_fuel] = fuel(FuelID, T);
R_univ = 8.31434e-3;      % kJ/mol-K
R_fuel = R_univ / (MW_fuel/1000); % kJ/kg-K (MW in kg/kmol → divide by 1000 to mol basis)
Cp_fuel = R_fuel * cp_fuel_nd;    % kJ/kg-K
h_fuel  = R_fuel * T * h_fuel_nd; % kJ/kg


if numel(varargin) == 2
    m_fa = varargin{1};
    m_res = varargin{2};

    m_fuel = m_fa / (AFR + 1);
    m_air  = m_fa - m_fuel;

    m_tot  = m_res + m_fa;
    h_fresh  = (m_air*h_air  + m_fuel*h_fuel ) / m_fa;                    % kJ/kg
    Cp_fresh = (m_air*Cp_air + m_fuel*Cp_fuel) / m_fa;                    % kJ/kg-K
    R_fresh  = (m_air*R_air  + m_fuel*R_fuel ) / m_fa;                    % kJ/kg-K
    h_mix   = (m_res*h_res   + m_fa*h_fresh ) / m_tot;               % kJ/kg
    Cp_mix  = (m_res*Cp_res  + m_fa*Cp_fresh) / m_tot;               % kJ/kg-K
    R_mix   = (m_res*R_res   + m_fa*R_fresh ) / m_tot;               % kJ/kg-K
    Cv_mix  = Cp_mix - R_mix;                                             % kJ/kg-K
else
    h_mix  = (AFR*h_air  + h_fuel ) / (AFR+1);                    % kJ/kg
    Cp_mix = (AFR*Cp_air + Cp_fuel) / (AFR+1);                    % kJ/kg-K
    R_mix  = (AFR*R_air  + R_fuel ) / (AFR+1);
    Cv_mix  = Cp_mix - R_mix;
end

end
