clear; clc; close all;

% Part Dimensions & Given Values
% ADD GIVEN VALUES HERE
Lv_max = 10;

b = 1;
b_dot = 1;

% Rise
% S = a4*(b^4) + a5*(b^5) + a6*(b^6) + a7*(b^7);
% V = 4*a4*(b^3)*b_dot + 5*a5*(b^4)*b_dot + 6*a6*(b^5)*b_dot + 7*a7*(b^6)*b_dot;
% A = 12*a4*(b^2)*(b_dot^2) + 20*a5*(b^3)*(b_dot^2) + 30*a6*(b^4)*(b_dot^2) + 42*a7*(b^5)*(b_dot^2);
% J = 24*a4*b*(b_dot^3) + 60*a5*(b^2)*(b_dot^3) + 120*a6*(b^3)*(b_dot^3) + 210*a7*(b^4)*(b_dot^3);

A = [b^4 b^5 b^6 b^7; 4*b^3 5*b^4 6*b^5 7*b^6; 12*b^2 20*b^3 30*b^4 42*b^5; 24*b 60*b^2 120*b^3 210*b^4];
B = [Lv_max; 0; 0; 0];
coefficients = A \ B;

a4 = coefficients(1);
a5 = coefficients(2);
a6 = coefficients(3);
a7 = coefficients(4);

x = 0:0.01:1;
for j = 1:101
    i = (j/100)-0.01;
    Lv(j) = a4*(i^4) + a5*(i^5) + a6*(i^6) + a7*(i^7);
end
plot(x,Lv);

%%
S_Eqn = a4*(x^4) + a5*(x^5) + a6*(x^6) + a7*(x^7);
V_Eqn = 4*a4*(x^3)*b_dot + 5*a5*(x^4)*b_dot + 6*a6*(x^5)*b_dot + 7*a7*(x^6)*b_dot;
A_Eqn = 12*a4*(x^2)*(b_dot^2) + 20*a5*(x^3)*(b_dot^2) + 30*a6*(x^4)*(b_dot^2) + 42*a7*(x^5)*(b_dot^2);
J_Eqn = 24*a4*x*(b_dot^3) + 60*a5*(x^2)*(b_dot^3) + 120*a6*(x^3)*(b_dot^3) + 210*a7*(x^4)*(b_dot^3);

% Return

% Dwell