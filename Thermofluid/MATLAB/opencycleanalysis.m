
% have some variable for start angle when valve opens and time that it
% takes to close
% have for intake and exhaust
% need some way to indicate how long valve is open and when
% leave a gap in the if statement that says when both valves are closed,
% call homogeneous function
%getting mdot in per second adjusted by degrees of crank angle, find what
%fraction of a second one degree of crank angle is
% function should set inputs as a variable and we'll figure out the inputs
% later, vars going to homogeneous are Tin, Pin, at end of intake need to
% have a temp and pressure val, also need the crank angle when switch over
% to homogeneous, overall outputs are heat and work
% run the valve function, will get bore and rpm from excel, will
% get intake and exhaust pressures and intake temperature from intake and
% exhaust functions.  wil also get wall temp from excel, for now just set
% it to 650k until we get a value
% remove line 78 in homogeneousV2, just use fixed value that it starts at
% will need to calc cp and cv inputs for the mixture using coolprop for air
% props and wiley fuel function for gas props, need these to feed into
% valve function
% use constant AFR of 14.7
% take cp and cv for the fuel and air and then take average based on either
% enthalpy, energy, or mass

function [heat,work] = opencycleanalysis(bore,stroke,rpm,T_intake,P_intake,P_exhaust)
clc

libLoc = "C\Users\Benjs\AppData\Roaming\CoolProp"

bore = bore;
rpm = rpm;
alpha_per_second = rpm/360;
stroke = stroke;
dim = [0.04128,45,0.037152,0.045408,0.01032]
fluid = 'Air'
lib = libLoc
AFR = 14.7;
%not sure what B is so i'm setting it to 1 for now
B = 1;

% not sure what phi is so i'm setting it as 0.9 for now
phi = 0.9;
ifuel = 2; %for gasoline

%starting at TDC with the intake valve closed
open = 0;
T_wall = 650;

% intake valve will open slightly at alpha = 1 degrees
alpha_intake_open_small = 1;
% intake valve will open partially at alpha = 10 degrees
alpha_intake_open_partial = 10;
% intake valve will open fully at alpha = 15 degrees
alpha_intake_open = 15;
% intake valve will close partially at alpha = 165 degrees
alpha_intake_closed_partial = 165;
% intake valve will close almost fully at alpha = 175 degrees
alpha_intake_closed_mostly = 175;
% intake valve will close at alpha = 180 degrees
alpha_intake_closed = 180;

% intake and exhaust valves both closed between alpha = 180 degrees and
% alpha = 539 degrees

% exhaust valve open slightly at alpha = 540 degrees
alpha_exhaust_open_small = 540;
% exhaust valve open partially at alpha = 545 degrees
alpha_exhaust_open_partial = 545;
% exhaust valve open fully at alpha = 555 degrees
alpha_exhaust_open = 555;
% exhaust valve closed partially at alpha = 705 degrees
alpha_exhaust_closed_partial = 705;
% exhaust valve closed mostly at alpha = 715 degrees
alpha_exhaust_closed_mostly = 715;
% exhaust valve closed fully at alpha = 720 degrees


% fluid props

for alpha = 1:720
    if alpha < alpha_open_partial
        % small intake valve opening at start of intake stroke
        %getting air props
        rho_inf(alpha) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(alpha) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(alpha) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        % setting valve open state
        open = 1;
        % getting the mdot for that alpha value
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(alpha) = mdot;
        % getting the fuel props
        [ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = ecp( T_intake, P_intake, phi, ifuel );
        cp_fuel(alpha) = Cp;
        cv_fuel(alpha) = v;
        % piston height based on the stroke minus the sin of 1/2 of the
        % crank angle multiplied by the stroke
        piston_height(alpha) = stroke - stroke*sind(alpha/2);
        % calculating the piston velocity
        deltay = piston_height(1);
        pistonvel(alpha) = deltay*alpha_per_second;
    elseif alpha >= alpha_open_partial & alpha < alpha_open
        % medium intake valve opening at start of intake stroke
        % getting air props
        rho_inf(alpha) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(alpha) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(alpha) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        % setting the valve open state
        open = 2;
        % getting the mdot for that alpha value
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(alpha) = mdot;
        % getting the fuel props
        [ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = ecp( T_intake, P_intake, phi, ifuel );
        cp_fuel(alpha) = Cp;
        cv_fuel(alpha) = v;
        % finding the piston height and velocity
        piston_height(alpha) = stroke - stroke*sind(alpha/2);
        pistonvel(alpha) = deltay*alpha_per_second;
    elseif alpha >= alpha_open & alpha < alpha_intake_closed_partial
        % fully open intake valve for most of intake stroke
        % getting the air props
        rho_inf(alpha) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(alpha) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(alpha) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        % setting the valve open state
        open = 3;
        % getting the mdot
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(alpha) = mdot;
        % getting the fuel props
        [ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = ecp( T_intake, P_intake, phi, ifuel );
        cp_fuel(alpha) = Cp;
        cv_fuel(alpha) = v;
        % finding the piston height and velocity
        piston_height(alpha) = stroke - stroke*sind(alpha/2);
        pistonvel(alpha) = deltay*alpha_per_second;
    elseif alpha >= alpha_intake_closed_partial & alpha < alpha_intake_closed_mostly
        % slightly closed intake valve at the end of the intake stroke
        rho_inf(alpha) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(alpha) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(alpha) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        open = 2;
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(alpha) = mdot;
        [ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = ecp( T_intake, P_intake, phi, ifuel );
        cp_fuel(alpha) = Cp;
        cv_fuel(alpha) = v;
        piston_height(alpha) = stroke - stroke*sind(alpha/2);
        pistonvel(alpha) = deltay*alpha_per_second;
    elseif alpha >= alpha_intake_closed_mostly & alpha < alpha_intake_closed
        % mostly closed intake valve at the end of the intake stroke
        rho_inf(alpha) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(alpha) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(alpha) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        open = 1;
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(alpha) = mdot;
        [ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = ecp( T_intake, P_intake, phi, ifuel );
        cp_fuel(alpha) = Cp;
        cv_fuel(alpha) = v;
        piston_height(alpha) = stroke - stroke*sind(alpha/2);
        pistonvel(alpha) = deltay*alpha_per_second;
    elseif alpha >= alpha_intake_closed & alpha < alpha_exhaust_open_small
        % both valves closed
        % homogeneous goes here
    elseif alpha >= alpha_exhaust_open_small & alpha < alpha_exhaust_open_partial
        % slightly open exhaust valve at the start of the exhaust stroke
        rho_inf(alpha) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(alpha) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(alpha) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        open = 1;
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(alpha) = mdot;
        [ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = ecp( T_intake, P_intake, phi, ifuel );
        cp_fuel(alpha) = Cp;
        cv_fuel(alpha) = v;
        piston_height(alpha) = stroke - stroke*abs(sind(alpha/2));
        deltay = piston_height(alpha_exhaust_open_small)
        pistonvel(alpha) = deltay*alpha_per_second;
    elseif alpha >= alpha_exhaust_open_partial & alpha < alpha_exhaust_open
        % mostly open exhaust valve at the start of the exhaust stroke
        rho_inf(alpha) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(alpha) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(alpha) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        open = 2;
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(alpha) = mdot;
        [ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = ecp( T_intake, P_intake, phi, ifuel );
        cp_fuel(alpha) = Cp;
        cv_fuel(alpha) = v;
        piston_height(alpha) = stroke - stroke*abs(sind(alpha/2));
        pistonvel(alpha) = deltay*alpha_per_second;
    elseif alpha >= alpha_exhaust_open & alpha < alpha_exhaust_closed_partial
        % fully open exhaust valve for most of the exhaust stroke
        rho_inf(alpha) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(alpha) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(alpha) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        open = 3;
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(alpha) = mdot;
        [ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = ecp( T_intake, P_intake, phi, ifuel );
        cp_fuel(alpha) = Cp;
        cv_fuel(alpha) = v;
        piston_height(alpha) = stroke - stroke*abs(sind(alpha/2));
        pistonvel(alpha) = deltay*alpha_per_second;
    elseif alpha >= alpha_exhaust_closed_partial & alpha < alpha_exhaust_closed_mostly
        % slightly closed exhaust valve for the end of the exhaust stroke
        rho_inf(alpha) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(alpha) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(alpha) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        open = 2;
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(alpha) = mdot;
        [ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = ecp( T_intake, P_intake, phi, ifuel );
        cp_fuel(alpha) = Cp;
        cv_fuel(alpha) = v;
        piston_height(alpha) = stroke - stroke*abs(sind(alpha/2));
        pistonvel(alpha) = deltay*alpha_per_second;
    elseif alpha >= alpha_exhaust_closed_mostly
        % mostly closed exhaust valve for the end of the exhaust stroke
        rho_inf(alpha) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(alpha) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(alpha) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        open = 1;
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(alpha) = mdot;
        [ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = ecp( T_intake, P_intake, phi, ifuel );
        cp_fuel(alpha) = Cp;
        cv_fuel(alpha) = v;
        piston_height(alpha) = stroke - stroke*abs(sind(alpha/2));
        pistonvel(alpha) = deltay*alpha_per_second;
    end

end

