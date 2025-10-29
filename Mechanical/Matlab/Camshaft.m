clear; clc; close all;

% Part Dimensions & Given Values
% ADD GIVEN VALUES HERE
Lv_max = 10;

% Rise
syms a4 a5 a6 a7
S = a4*(beta^4) + a5*(beta^5) + a6*(beta^6) + a7*(beta^7);
V = 4*a4*(beta^3)*beta_dot + 5*a5*(beta^4)*beta_dot + 6*a6*(beta^5)*beta_dot + 7*a7*(beta^6)*beta_dot;
A = 12*a4*(beta^2)*(beta_dot^2) + 20*a5*(beta^3)*(beta_dot^2) + 30*a6*(beta^4)*(beta_dot^2) + 42*a7*(beta^5)*(beta_dot^2);
J = 24*a4*beta*(beta_dot^3) + 60*a5*(beta^2)*(beta_dot^3) + 120*a6*(beta^3)*(beta_dot^3) + 210*a7*(beta^4)*(beta_dot^3);

[a4,a5,a6,a7] = solve(S(beta=1) == Lv_max, V(beta=1, beta_dot="IDK") == 0, A(beta=1, beta_dot="IDK") == 0, J(beta=1, beta_dot="IDK") == 0);

% Return

% Dwell