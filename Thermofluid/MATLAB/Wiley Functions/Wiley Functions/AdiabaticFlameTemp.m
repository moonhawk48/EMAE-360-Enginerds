
% Computes const pressure adiabatic flame temperature T2 
% from first law: dh = q
% Inputs:
  T1 = 600;% initial temperature (K)
  P1 = 100; % initial pressure (kPa)
  PHI = 1.0; %  equivalence ratio
  f = .1; %  residual fraction
  ifuel=2; % 1=Methane, 2=Gasoline, 3=Diesel, 4=Methanol, 5=Nitromethane
% use FARG for initial calc of H1
[~, H1, ~, ~, ~, ~, CP, ~, ~, ~] = farg( T1, P1, PHI, f, ifuel );
fprintf(' Initial Enthalpy H1 =  %7.1f Initial CP = %7.3f \n',H1,CP);
P2 = P1;
T2 = 2000; %initial guess of flame temp
MAXITS = 50;
TOL = 0.00001;

for i=1:MAXITS,
    [~,~, H2, ~, ~, ~, ~, CP, ~, ~, ~] = ecp( T2, P2, PHI, ifuel );
    %fprintf(' Iterated Enthalpy H2 =  %7.1f Iterated CP = %7.3f \n',H2,CP);
    DELT2 = (H1-H2)/CP;
    T2 = T2 + DELT2;
    fprintf(' Iterated Adiabatic Flame Temp (K) =  %8.1f  \n',T2);
    if ( abs(DELT2)/T2 < TOL )
        break;
    end
end
fprintf(' Final Adiabatic Flame Temp (K) =  %8.1f  \n',T2);