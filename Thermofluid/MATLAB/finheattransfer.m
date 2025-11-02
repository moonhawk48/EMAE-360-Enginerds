

function [N] = fincool(L,hair,Tair,Qin,Tb,tfin,Ab,k,P,Rb)
theta = Tb-Tair;
Ac = pi*(L^2) - pi*(Rb^2);
m = sqrt((hair*P)/(k*Ac));
heat_flow_per_temp_gradient = Qin/theta;
N_with_extras = heat_flow_per_temp_gradient - hair*Ab;
a = tanh(m*L);
N = N_with_extras/(a*sqrt(hair*P*k*Ac));

end