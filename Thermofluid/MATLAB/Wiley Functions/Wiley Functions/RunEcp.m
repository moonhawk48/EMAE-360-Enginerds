clear;
phi = 0.9;  % enter equivalence ratio input
T = 800; % enter temperature (K) input
P = 100;  % enter pressure (kPa) input
fuel_id =1;  
%   fuel_id - 1=Methane, 2=Gasoline, 3=Diesel, 4=Methanol, 5=Nitromethane
% call ecp function
[ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = ecp( T, P, phi, fuel_id );

%echo input 
fprintf(' \n Equilibrium Combustion Solver \n' );
fprintf(' Pressure (kPa) = \t \t %6.1f \n', P );
fprintf(' Temperature (K) = \t \t %6.1f \n', T); ...
fprintf(' Fuel Air Equivalence ratio = \t% 4.2f \n ', phi);
 
%print output mole fractions and properties
fprintf(' \n Mole Fractions \n' );
fprintf(' CO2 = \t %1.4e \n', Y(1) );
fprintf(' H2O = \t %1.4e \n', Y(2) );
fprintf(' N2 = \t %1.4e \n', Y(3) );
fprintf(' O2 = \t %1.4e \n', Y(4) );
fprintf(' CO = \t %1.4e \n', Y(5) );
fprintf(' H2 = \t %1.4e \n', Y(6) );
fprintf(' H = \t %1.4e \n', Y(7) );
fprintf(' O = \t %1.4e \n', Y(8) );
fprintf(' OH = \t %1.4e \n', Y(9) );
fprintf(' NO = \t %1.4e \n', Y(10) );

fprintf(' \n Mixture Properties \n' );
fprintf(' h(kJ/kg) = \t %6.1f \n', h );
fprintf(' u(kJ/kg) = \t %6.1f \n', u );
fprintf(' s (kJ/Kg K) = \t %6.4f \n', s );
fprintf(' v (m3/kg) = \t %6.4f \n', v );
fprintf(' cp (kJ/Kg K) =\t %6.3f \n', Cp );
fprintf(' Molecular Mass = %5.2f \n', MW );
fprintf(' dvdt = %1.4e \n', dvdT );
fprintf(' dvdp = %1.4e \n', dvdP );