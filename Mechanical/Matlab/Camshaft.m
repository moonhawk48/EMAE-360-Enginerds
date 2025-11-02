clear; clc; close all;

% Part Dimensions & Given Values
% ADD MORE PART DIMENSIONS HERE
Lv_max = 16.9; % Max Valve Lift [mm]

b = 1;
b_dot = 1;

% Solving for Coefficients of 4-5-6-7 Polynomial
A = [b^4 b^5 b^6 b^7; 4*b^3 5*b^4 6*b^5 7*b^6; 12*b^2 20*b^3 30*b^4 42*b^5; 24*b 60*b^2 120*b^3 210*b^4];
B = [Lv_max; 0; 0; 0];
coefficients = A \ B;

a4 = coefficients(1);
a5 = coefficients(2);
a6 = coefficients(3);
a7 = coefficients(4);

% --------- Plotting SVAJ Curves ---------
% NOT SURE ABOUT THE AXES EXACTLY SO NEED TO CHECK THEM

% Displacement
x = 0:0.01:1;
for j = 1:101
    i = (j/100)-0.01;
    Lv(j) = a4*(i^4) + a5*(i^5) + a6*(i^6) + a7*(i^7);
end
subplot(4,1,1);
sgtitle("SVAJ Curves");
plot(x,Lv);
xlabel("Crank Angle [θ]");
ylabel("S - Displacement [mm]");

%  Velocity
for j = 1:101
    i = (j/100)-0.01;
    Lv_dot(j) = 4*a4*(i^3)*b_dot + 5*a5*(i^4)*b_dot + 6*a6*(i^5)*b_dot + 7*a7*(i^6)*b_dot;
end
subplot(4,1,2);
plot(x,Lv_dot);
xlabel("Crank Angle [θ]");
ylabel("V - Velocity [mm/s]");

%  Acceleration
for j = 1:101
    i = (j/100)-0.01;
    Lv_ddot(j) = 12*a4*(i^2)*(b_dot^2) + 20*a5*(i^3)*(b_dot^2) + 30*a6*(i^4)*(b_dot^2) + 42*a7*(i^5)*(b_dot^2);
end
subplot(4,1,3);
plot(x,Lv_ddot);
xlabel("Crank Angle [θ]");
ylabel("A - Acceleration [mm/s^2]");

%  Jerk
for j = 1:101
    i = (j/100)-0.01;
    Lv_tdot(j) = 24*a4*i*(b_dot^3) + 60*a5*(i^2)*(b_dot^3) + 120*a6*(i^3)*(b_dot^3) + 210*a7*(i^4)*(b_dot^3);
end
subplot(4,1,4);
plot(x,Lv_tdot);
xlabel("Crank Angle [θ]");
ylabel("J - Jerk [mm/s^3]");

%% Trying Reverse Plots for Return Motion
% Displacement
x = 1:0.01:2;
for j = 1:101
    i = (j/100)-0.01;
    Lv(j) = a4*(i^4) + a5*(i^5) + a6*(i^6) + a7*(i^7);
end
figure(2);
subplot(4,1,1);
sgtitle("SVAJ Curves");
plot(x,fliplr(Lv));
xlabel("Crank Angle [θ]");
ylabel("S - Displacement [mm]");

%  Velocity
for j = 1:101
    i = (j/100)-0.01;
    Lv_dot(j) = 4*a4*(i^3)*b_dot + 5*a5*(i^4)*b_dot + 6*a6*(i^5)*b_dot + 7*a7*(i^6)*b_dot;
end
subplot(4,1,2);
plot(x,fliplr(Lv_dot));
xlabel("Crank Angle [θ]");
ylabel("V - Velocity [mm/s]");

%  Acceleration
for j = 1:101
    i = (j/100)-0.01;
    Lv_ddot(j) = 12*a4*(i^2)*(b_dot^2) + 20*a5*(i^3)*(b_dot^2) + 30*a6*(i^4)*(b_dot^2) + 42*a7*(i^5)*(b_dot^2);
end
subplot(4,1,3);
plot(x,fliplr(Lv_ddot));
xlabel("Crank Angle [θ]");
ylabel("A - Acceleration [mm/s^2]");

%  Jerk
for j = 1:101
    i = (j/100)-0.01;
    Lv_tdot(j) = 24*a4*i*(b_dot^3) + 60*a5*(i^2)*(b_dot^3) + 120*a6*(i^3)*(b_dot^3) + 210*a7*(i^4)*(b_dot^3);
end
subplot(4,1,4);
plot(x,fliplr(Lv_tdot));
xlabel("Crank Angle [θ]");
ylabel("J - Jerk [mm/s^3]");


%%
% NEXT, NEED TO PLOT THE SHAPE OF THE CAM PROFILE USING X & Y MAPPING FROM
% BACHMANN TEXTBOOK