
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

function [mdottotal,Etotal,Qtotal,Wtotal,Htotal,T_exhaust_pipe_total,P_exhaust_pipe_total] = opencycleanalysis(bore,stroke,rpm,T_intake,P_intake,P_exhaust)
clc

libLoc = "C\Users\Benjs\AppData\Roaming\CoolProp"

bore = bore;
rpm = rpm;
alpha_per_second = (rpm/360)/60;
seconds_per_alpha = 1/alpha_per_second;
stroke = stroke;
dim = [0.04128,45,0.037152,0.045408,0.01032]
fluid = 'Air'
lib = libLoc
AFR = 14.7;
%not sure what B is so i'm setting it to 1 for now
B = 1;


id = 2; %for gasoline
phi = 1/AFR;
fuel_id = 2;

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


cylinder_area = pi*(bore/2)^2;
% fluid props

rho_inf = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
cp_air = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
cv_air = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
[ angle, beta, gamma, delta, h, s, cp, mw, Fs, q ] = fuel( id, T_intake);
cp_fuel = cp;
cv_fuel = cp/gamma;
cp_total = cp_fuel+cp_air;
cv_total = cv_fuel+cv_air;
gam_total = (cp_fuel+AFR*cp_air)/(cv_fuel+AFR*cv_air);

q_init = q;


Tend = 1350;
Pend = 101325*3; % in Pa, not sure what units matlab stuff uses
f = 0.1;
[Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = farg( Tend, Pend, phi, f, fuel_id )
h_comb = h;
Eend = u*Tend;

delP = Pend-P_intake;
hq = woschni(stage,B,Pend,Tend,MPS,T_intake,V,Vo,delP,Pend)

Q = hq*cylinder_area*(Tend-T_wall);

Hend = mw*cp_total*Tend;



%ballpark Tend and Pend Tend is somewhere between 1200 and 1500k, Pend is
%near atmosphere, maybe 3 atm
% use farg wiley function to get specific enthalpy and energy at this state

% to get q, use the woschni function, just use the intake and
% exhaust cases to get h
% use this function: Q = h*Cylinder Area(Tgas-Twall) to get Q
% use internal energy to get the temperature of the gas, wall temp will
% later get passed into it
% Qinto air will always be negative in this section
% E_new = E_old - Q + W - H
% W=P*V=P*A*deltah
% Winto air will always be positive
% H = m*(cp)mixture*T


for angle = 1:720
    if angle < alpha_open_partial
        % small intake valve opening at start of intake stroke
        %getting air props
        rho_inf(angle) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(angle) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(angle) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        % setting valve open state
        open = 1;
        % getting the mdot for that alpha value
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(angle) = mdot;
        % getting the new temperature
        if angle ==1
            q(angle) = q_init;
            cp_fuel(angle) = cp_fuel;
        else
            q(angle) = q(angle-1);
            cp_fuel(angle) = cp_fuel(angle-1);
        end
        deltaT(angle) = q(angle)/(mdotalpha(angle)*cp_fuel(angle));
        if angle ==1
            T(angle) = Tend;
        else
            T(angle) = T(angle-1) + deltaT(angle);
        end

        % getting the fuel props
        [ alpha, beta, gamma, delta, h, s, cp, mw, Fs, q ] = fuel( id, T(angle));
        mw(angle) = mw;
        cp_fuel(angle) = cp;
        cv_fuel(angle) = cp/gamma;
        cp_total(angle) = cp_fuel(angle)+cp_air(angle);
        cv_total(angle) = cv_fuel(angle)+cv_air(angle);
        gam_total(angle) = (cp_fuel(angle)+AFR*cp_air(angle))/(cv_fuel(angle)+AFR*cv_air(angle));
        P_exhaust_pipe(angle) = 0;
        T_exhaust_pipe(angle) = 0;
        % finding the internal energy
        P(angle) = getFluidProperty(libLoc,'P','CVMASS',cv_air(angle),'T',T(angle));
        [Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = farg( T(angle), T(angle), phi, f, fuel_id );
        h_tot(angle) = h;
        E(angle) = u*T(angle);
        delP(angle) = P(angle)-P_intake;
        hq(angle) = woschni(stage,B,P(angle),T(angle),MPS,T_intake,V,Vo,delP(angle),P_intake);
        Q(angle) = hq(angle)*cylinder_area*(T(angle)-T_wall);
        H(angle) = mw(angle)*cp_total(angle)*T(angle);
        
        % piston height based on the stroke minus the sin of 1/2 of the
        % crank angle multiplied by the stroke
        piston_height(angle) = stroke - stroke*sind(angle/2);
        % calculating the piston velocity
        deltay(angle) = stroke*sind(angle/2);
        pistonvel(angle) = deltay*alpha_per_second;
        % getting the work
        w(angle) = P_intake*cylinder_area*abs(deltay(angle));
    elseif angle >= alpha_open_partial & angle < alpha_open
        % medium intake valve opening at start of intake stroke
        % getting air props
        rho_inf(angle) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(angle) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(angle) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        % setting the valve open state
        open = 2;
        % getting the mdot for that alpha value
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(angle) = mdot;
        % getting the new temperature
        if angle ==1
            q(angle) = q_init;
            cp_fuel(angle) = cp_fuel;
        else
            q(angle) = q(angle-1);
            cp_fuel(angle) = cp_fuel(angle-1);
        end
        deltaT(angle) = q(angle)/(mdotalpha(angle)*cp_fuel(angle));
        if angle ==1
            T(angle) = Tend;
        else
            T(angle) = T(angle-1) + deltaT(angle);
        end
        % getting the fuel props
        [ alpha, beta, gamma, delta, h, s, cp, mw, Fs, q ] = fuel( id, T(angle));
        mw(angle) = mw;
        cp_fuel(angle) = cp;
        cv_fuel(angle) = cp/gamma;
        cp_total(angle) = cp_fuel(angle)+cp_air(angle);
        cv_total(angle) = cv_fuel(angle)+cv_air(angle);
        gam_total(angle) = (cp_fuel(angle)+AFR*cp_air(angle))/(cv_fuel(angle)+AFR*cv_air(angle));
        P_exhaust_pipe(angle) = 0;
        T_exhaust_pipe(angle) = 0;
        % finding the internal energy
        P(angle) = getFluidProperty(libLoc,'P','CVMASS',cv_air(angle),'T',T(angle));
        [Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = farg( T(angle), T(angle), phi, f, fuel_id );
        h_tot(angle) = h;
        E(angle) = u*T(angle);
        delP(angle) = P(angle)-P_intake;
        hq(angle) = woschni(stage,B,P(angle),T(angle),MPS,T_intake,V,Vo,delP(angle),P_intake);
        Q(angle) = hq(angle)*cylinder_area*(T(angle)-T_wall);
        H(angle) = mw(angle)*cp_total(angle)*T(angle);
        % finding the piston height and velocity
        piston_height(angle) = stroke - stroke*sind(angle/2);
        deltay(angle) = stroke*sind(angle/2);
        pistonvel(angle) = deltay*alpha_per_second;
        % finding the work
        w(angle) = P_intake*cylinder_area*abs(deltay(angle));
    elseif angle >= alpha_open & angle < alpha_intake_closed_partial
        % fully open intake valve for most of intake stroke
        % getting the air props
        rho_inf(angle) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(angle) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(angle) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        % setting the valve open state
        open = 3;
        % getting the mdot
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(angle) = mdot;
        % getting the new temperature
        if angle ==1
            q(angle) = q_init;
            cp_fuel(angle) = cp_fuel;
        else
            q(angle) = q(angle-1);
            cp_fuel(angle) = cp_fuel(angle-1);
        end
        deltaT(angle) = q(angle)/(mdotalpha(angle)*cp_fuel(angle));
        if angle ==1
            T(angle) = Tend;
        else
            T(angle) = T(angle-1) + deltaT(angle);
        end
        % getting the fuel props
        [ alpha, beta, gamma, delta, h, s, cp, mw, Fs, q ] = fuel( id, T(angle));
        mw(angle) = mw;
        cp_fuel(angle) = cp;
        cv_fuel(angle) = cp/gamma;
        cp_total(angle) = cp_fuel(angle)+cp_air(angle);
        cv_total(angle) = cv_fuel(angle)+cv_air(angle);
        gam_total(angle) = (cp_fuel(angle)+AFR*cp_air(angle))/(cv_fuel(angle)+AFR*cv_air(angle));
        P_exhaust_pipe(angle) = 0;
        T_exhaust_pipe(angle) = 0;
        % finding the internal energy
        P(angle) = getFluidProperty(libLoc,'P','CVMASS',cv_air(angle),'T',T(angle));
        [Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = farg( T(angle), T(angle), phi, f, fuel_id );
        h_tot(angle) = h;
        E(angle) = u*T(angle);
        delP(angle) = P(angle)-P_intake;
        hq(angle) = woschni(stage,B,P(angle),T(angle),MPS,T_intake,V,Vo,delP(angle),P_intake);
        Q(angle) = hq(angle)*cylinder_area*(T(angle)-T_wall);
        H(angle) = mw(angle)*cp_total(angle)*T(angle);
        % finding the piston height and velocity
        piston_height(angle) = stroke - stroke*sind(angle/2);
        deltay(angle) = stroke*sind(angle/2);
        pistonvel(angle) = deltay*alpha_per_second;
        % finding the work
        w(angle) = P_intake*cylinder_area*abs(deltay(angle));
    elseif angle >= alpha_intake_closed_partial & angle < alpha_intake_closed_mostly
        % slightly closed intake valve at the end of the intake stroke
        % getting the air props
        rho_inf(angle) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(angle) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(angle) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        % setting the valve open state
        open = 2;
        % getting the mdot
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(angle) = mdot;
        % getting the new temperature
        if angle ==1
            q(angle) = q_init;
            cp_fuel(angle) = cp_fuel;
        else
            q(angle) = q(angle-1);
            cp_fuel(angle) = cp_fuel(angle-1);
        end
        deltaT(angle) = q(angle)/(mdotalpha(angle)*cp_fuel(angle));
        if angle ==1
            T(angle) = Tend;
        else
            T(angle) = T(angle-1) + deltaT(angle);
        end
        % getting the fuel props
        [ alpha, beta, gamma, delta, h, s, cp, mw, Fs, q ] = fuel( id, T(angle));
        mw(angle) = mw;
        cp_fuel(angle) = cp;
        cv_fuel(angle) = cp/gamma;
        cp_total(angle) = cp_fuel(angle)+cp_air(angle);
        cv_total(angle) = cv_fuel(angle)+cv_air(angle);
        gam_total(angle) = (cp_fuel(angle)+AFR*cp_air(angle))/(cv_fuel(angle)+AFR*cv_air(angle));
        P_exhaust_pipe(angle) = 0;
        T_exhaust_pipe(angle) = 0;
        % finding the internal energy
        P(angle) = getFluidProperty(libLoc,'P','CVMASS',cv_air(angle),'T',T(angle));
        [Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = farg( T(angle), T(angle), phi, f, fuel_id );
        h_tot(angle) = h;
        E(angle) = u*T(angle);
        delP(angle) = P(angle)-P_intake;
        hq(angle) = woschni(stage,B,P(angle),T(angle),MPS,T_intake,V,Vo,delP(angle),P_intake);
        Q(angle) = hq(angle)*cylinder_area*(T(angle)-T_wall);
        H(angle) = mw(angle)*cp_total(angle)*T(angle);
        % finding the piston height and velocity
        piston_height(angle) = stroke - stroke*sind(angle/2);
        deltay(angle) = stroke*sind(angle/2);
        pistonvel(angle) = deltay*alpha_per_second;
        % finding the work
        w(angle) = P_intake*cylinder_area*abs(deltay(angle));
    elseif angle >= alpha_intake_closed_mostly & angle < alpha_intake_closed
        % mostly closed intake valve at the end of the intake stroke
        % getting the air props
        rho_inf(angle) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(angle) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(angle) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        % setting the valve open state
        open = 1;
        % getting the mdot
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(angle) = mdot;
        % getting the new temperature
        if angle ==1
            q(angle) = q_init;
            cp_fuel(angle) = cp_fuel;
        else
            q(angle) = q(angle-1);
            cp_fuel(angle) = cp_fuel(angle-1);
        end
        deltaT(angle) = q(angle)/(mdotalpha(angle)*cp_fuel(angle));
        if angle ==1
            T(angle) = Tend;
        else
            T(angle) = T(angle-1) + deltaT(angle);
        end
        % getting the fuel props
        [ alpha, beta, gamma, delta, h, s, cp, mw, Fs, q ] = fuel( id, T(angle));
        mw(angle) = mw;
        cp_fuel(angle) = cp;
        cv_fuel(angle) = cp/gamma;
        cp_total(angle) = cp_fuel(angle)+cp_air(angle);
        cv_total(angle) = cv_fuel(angle)+cv_air(angle);
        gam_total(angle) = (cp_fuel(angle)+AFR*cp_air(angle))/(cv_fuel(angle)+AFR*cv_air(angle));
        P_exhaust_pipe(angle) = 0;
        T_exhaust_pipe(angle) = 0;
        % finding the internal energy
        P(angle) = getFluidProperty(libLoc,'P','CVMASS',cv_air(angle),'T',T(angle));
        [Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = farg( T(angle), T(angle), phi, f, fuel_id );
        h_tot(angle) = h;
        E(angle) = u*T(angle);
        delP(angle) = P(angle)-P_intake;
        hq(angle) = woschni(stage,B,P(angle),T(angle),MPS,T_intake,V,Vo,delP(angle),P_intake);
        Q(angle) = hq(angle)*cylinder_area*(T(angle)-T_wall);
        H(angle) = mw(angle)*cp_total(angle)*T(angle);
        % finding the piston height and velocity
        piston_height(angle) = stroke - stroke*sind(angle/2);
        deltay(angle) = stroke*sind(angle/2);
        pistonvel(angle) = deltay*alpha_per_second;
        % finding the work
        w(angle) = P_intake*cylinder_area*abs(deltay(angle));
    elseif angle >= alpha_intake_closed & angle < alpha_exhaust_open_small
        % both valves closed
        % homogeneous goes here
        % include T and P at end of intake, also send wall temp and total
        % charge mass: sum of the mass taken in at each crank angle: do
        % this before you call the closed system function
        compmasstot = sum(mdot);
        varagin = [T(angle-1),P(angle-1),RPM,T_wall,compmasstot];
        [W_total,Q_total,T_out,P_out] = Homogeneous(varargin);
        W_total(angle) = W_total;
        w(angle) = W_total(angle);
        Q_total(angle) = Q_total;
        Q(angle) = Q_total(angle);
        T_out(angle) = T_out;
        T(angle) = T_out(angle);
        T_exhaust_pipe(angle) = T(angle);
        P_out(angle) = P_out;
        P(angle) = P_out(angle);
        P_exhaust_pipe(angle) = P(angle);
    elseif angle >= alpha_exhaust_open_small & angle < alpha_exhaust_open_partial
        % slightly open exhaust valve at the start of the exhaust stroke
        % getting the air props
        rho_inf(angle) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(angle) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(angle) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        % setting the valve open state
        open = 1;
        % getting the mdot
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(angle) = mdot;
        % getting the new temperature
        if angle ==1
            q(angle) = q_init;
            cp_fuel(angle) = cp_fuel;
        else
            q(angle) = q(angle-1);
            cp_fuel(angle) = cp_fuel(angle-1);
        end
        deltaT(angle) = q(angle)/(mdotalpha(angle)*cp_fuel(angle));
        if angle ==1
            T(angle) = Tend;
        else
            T(angle) = T(angle-1) + deltaT(angle);
        end
        % getting the air props
        [ alpha, beta, gamma, delta, h, s, cp, mw, Fs, q ] = fuel( id, T(angle));
        mw(angle) = mw;
        cp_fuel(angle) = cp;
        cv_fuel(angle) = cp/gamma;
        cp_total(angle) = cp_fuel(angle)+cp_air(angle);
        cv_total(angle) = cv_fuel(angle)+cv_air(angle);
        gam_total(angle) = (cp_fuel(angle)+AFR*cp_air(angle))/(cv_fuel(angle)+AFR*cv_air(angle));
        % finding the internal energy
        P(angle) = getFluidProperty(libLoc,'P','CVMASS',cv_air(angle),'T',T(angle));
        P_exhaust_pipe(angle) = P(angle);
        T_exhaust_pipe(angle) = T(angle);
        [Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = farg( T(angle), T(angle), phi, f, fuel_id );
        h_tot(angle) = h;
        E(angle) = u*T(angle);
        delP(angle) = P(angle)-P_intake;
        hq(angle) = woschni(stage,B,P(angle),T(angle),MPS,T_intake,V,Vo,delP(angle),P_intake);
        Q(angle) = hq(angle)*cylinder_area*(T(angle)-T_wall);
        H(angle) = mw(angle)*cp_total(angle)*T(angle);
        % finding the piston height and velocity
        piston_height(angle) = stroke - stroke*abs(sind(angle/2));
        deltay = piston_height(alpha_exhaust_open_small);
        deltay(angle) = stroke*sind(angle/2);
        pistonvel(angle) = deltay*alpha_per_second;
        % finding the work
        w(angle) = P_intake*cylinder_area*abs(deltay(angle));
    elseif angle >= alpha_exhaust_open_partial & angle < alpha_exhaust_open
        % mostly open exhaust valve at the start of the exhaust stroke
        % getting the air props
        rho_inf(angle) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(angle) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(angle) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        % setting the valve open state
        open = 2;
        %getting the mdot
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(angle) = mdot;
        % getting the new temperature
        if angle ==1
            q(angle) = q_init;
            cp_fuel(angle) = cp_fuel;
        else
            q(angle) = q(angle-1);
            cp_fuel(angle) = cp_fuel(angle-1);
        end
        deltaT(angle) = q(angle)/(mdotalpha(angle)*cp_fuel(angle));
        if angle ==1
            T(angle) = Tend;
        else
            T(angle) = T(angle-1) + deltaT(angle);
        end
        % getting the fuel props
        [ alpha, beta, gamma, delta, h, s, cp, mw, Fs, q ] = fuel( id, T(angle));
        mw(angle) = mw;
        cp_fuel(angle) = cp;
        cv_fuel(angle) = cp/gamma;
        cp_total(angle) = cp_fuel(angle)+cp_air(angle);
        cv_total(angle) = cv_fuel(angle)+cv_air(angle);
        gam_total(angle) = (cp_fuel(angle)+AFR*cp_air(angle))/(cv_fuel(angle)+AFR*cv_air(angle));
        % finding the internal energy
        P(angle) = getFluidProperty(libLoc,'P','CVMASS',cv_air(angle),'T',T(angle));
        P_exhaust_pipe(angle) = P(angle);
        T_exhaust_pipe(angle) = T(angle);
        [Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = farg( T(angle), T(angle), phi, f, fuel_id );
        h_tot(angle) = h;
        E(angle) = u*T(angle);
        delP(angle) = P(angle)-P_intake;
        hq(angle) = woschni(stage,B,P(angle),T(angle),MPS,T_intake,V,Vo,delP(angle),P_intake);
        Q(angle) = hq(angle)*cylinder_area*(T(angle)-T_wall);
        H(angle) = mw(angle)*cp_total(angle)*T(angle);
        % finding the piston height and velocity
        piston_height(angle) = stroke - stroke*abs(sind(angle/2));
        deltay(angle) = stroke*sind(angle/2);
        pistonvel(angle) = deltay*alpha_per_second;
        % finding the work
        w(angle) = P_intake*cylinder_area*abs(deltay(angle));
    elseif angle >= alpha_exhaust_open & angle < alpha_exhaust_closed_partial
        % fully open exhaust valve for most of the exhaust stroke
        % getting the air props
        rho_inf(angle) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(angle) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(angle) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        % setting the valve open state
        open = 3;
        % getting the mdot
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(angle) = mdot;
        % getting the new temperature
        if angle ==1
            q(angle) = q_init;
            cp_fuel(angle) = cp_fuel;
        else
            q(angle) = q(angle-1);
            cp_fuel(angle) = cp_fuel(angle-1);
        end
        deltaT(angle) = q(angle)/(mdotalpha(angle)*cp_fuel(angle));
        if angle ==1
            T(angle) = Tend;
        else
            T(angle) = T(angle-1) + deltaT(angle);
        end
        % getting the fuel props
        [ alpha, beta, gamma, delta, h, s, cp, mw, Fs, q ] = fuel( id, T(angle));
        mw(angle) = mw;
        cp_fuel(angle) = cp;
        cv_fuel(angle) = cp/gamma;
        cp_total(angle) = cp_fuel(angle)+cp_air(angle);
        cv_total(angle) = cv_fuel(angle)+cv_air(angle);
        gam_total(angle) = (cp_fuel(angle)+AFR*cp_air(angle))/(cv_fuel(angle)+AFR*cv_air(angle));
        % finding the internal energy
        P(angle) = getFluidProperty(libLoc,'P','CVMASS',cv_air(angle),'T',T(angle));
        P_exhaust_pipe(angle) = P(angle);
        T_exhaust_pipe(angle) = T(angle);
        [Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = farg( T(angle), T(angle), phi, f, fuel_id );
        h_tot(angle) = h;
        E(angle) = u*T(angle);
        delP(angle) = P(angle)-P_intake;
        hq(angle) = woschni(stage,B,P(angle),T(angle),MPS,T_intake,V,Vo,delP(angle),P_intake);
        Q(angle) = hq(angle)*cylinder_area*(T(angle)-T_wall);
        H(angle) = mw(angle)*cp_total(angle)*T(angle);
        % finding the piston height and velocity
        piston_height(angle) = stroke - stroke*abs(sind(angle/2));
        deltay(angle) = stroke*sind(angle/2);
        pistonvel(angle) = deltay*alpha_per_second;
        % finding the work
        w(angle) = P_intake*cylinder_area*abs(deltay(angle));
    elseif angle >= alpha_exhaust_closed_partial & angle < alpha_exhaust_closed_mostly
        % slightly closed exhaust valve for the end of the exhaust stroke
        % getting the air props
        rho_inf(angle) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(angle) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(angle) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        % setting the valve open state
        open = 2;
        % getting the mdot
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(angle) = mdot;
        % getting the new temperature
        if angle ==1
            q(angle) = q_init;
            cp_fuel(angle) = cp_fuel;
        else
            q(angle) = q(angle-1);
            cp_fuel(angle) = cp_fuel(angle-1);
        end
        deltaT(angle) = q(angle)/(mdotalpha(angle)*cp_fuel(angle));
        if angle ==1
            T(angle) = Tend;
        else
            T(angle) = T(angle-1) + deltaT(angle);
        end
        % getting the fuel props
        [ alpha, beta, gamma, delta, h, s, cp, mw, Fs, q ] = fuel( id, T(angle));
        mw(angle) = mw;
        cp_fuel(angle) = cp;
        cv_fuel(angle) = cp/gamma;
        cp_total(angle) = cp_fuel(angle)+cp_air(angle);
        cv_total(angle) = cv_fuel(angle)+cv_air(angle);
        gam_total(angle) = (cp_fuel(angle)+AFR*cp_air(angle))/(cv_fuel(angle)+AFR*cv_air(angle));
        % finding the internal energy
        P(angle) = getFluidProperty(libLoc,'P','CVMASS',cv_air(angle),'T',T(angle));
        P_exhaust_pipe(angle) = P(angle);
        T_exhaust_pipe(angle) = T(angle);
        [Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = farg( T(angle), T(angle), phi, f, fuel_id );
        h_tot(angle) = h;
        E(angle) = u*T(angle);
        delP(angle) = P(angle)-P_intake;
        hq(angle) = woschni(stage,B,P(angle),T(angle),MPS,T_intake,V,Vo,delP(angle),P_intake);
        Q(angle) = hq(angle)*cylinder_area*(T(angle)-T_wall);
        H(angle) = mw(angle)*cp_total(angle)*T(angle);
        % finding the piston height and velocity
        piston_height(angle) = stroke - stroke*abs(sind(angle/2));
        deltay(angle) = stroke*sind(angle/2);
        pistonvel(angle) = deltay*alpha_per_second;
        % finding the work
        w(angle) = P_intake*cylinder_area*abs(deltay(angle));
    elseif angle >= alpha_exhaust_closed_mostly
        % mostly closed exhaust valve for the end of the exhaust stroke
        % getting the air props
        rho_inf(angle) = getFluidProperty(libLoc,'D','T',T_intake,'P',P_intake,fluid);
        cp_air(angle) = getFluidProperty(libLoc,'CPMASS','T',T_intake,'P',P_intake,fluid);
        cv_air(angle) = getFluidProperty(libLoc,'CVMASS','T',T_intake,'P',P_intake,fluid);
        % setting the valve open state
        open = 1;
        % getting the mdot
        mdot = valve(open,B,dim,T_intake,P_intake,P_exhaust,fluid,lib);
        mdotalpha(angle) = mdot;
        % getting the new temperature
        if angle ==1
            q(angle) = q_init;
            cp_fuel(angle) = cp_fuel;
        else
            q(angle) = q(angle-1);
            cp_fuel(angle) = cp_fuel(angle-1);
        end
        deltaT(angle) = q(angle)/(mdotalpha(angle)*cp_fuel(angle));
        if angle ==1
            T(angle) = Tend;
        else
            T(angle) = T(angle-1) + deltaT(angle);
        end
        % getting the fuel props
        [ alpha, beta, gamma, delta, h, s, cp, mw, Fs, q ] = fuel( id, T(angle));
        mw(angle) = mw;
        cp_fuel(angle) = cp;
        cv_fuel(angle) = cp/gamma;
        cp_total(angle) = cp_fuel(angle)+cp_air(angle);
        cv_total(angle) = cv_fuel(angle)+cv_air(angle);
        gam_total(angle) = (cp_fuel(angle)+AFR*cp_air(angle))/(cv_fuel(angle)+AFR*cv_air(angle));
        % finding the internal energy
        P(angle) = getFluidProperty(libLoc,'P','CVMASS',cv_air(angle),'T',T(angle));
        P_exhaust_pipe(angle) = P(angle);
        T_exhaust_pipe(angle) = T(angle);
        [Y, h, u, s, v, R, Cp, MW, dvdT, dvdP] = farg( T(angle), T(angle), phi, f, fuel_id );
        h_tot(angle) = h;
        E(angle) = u*T(angle);
        delP(angle) = P(angle)-P_intake;
        hq(angle) = woschni(stage,B,P(angle),T(angle),MPS,T_intake,V,Vo,delP(angle),P_intake);
        Q(angle) = hq(angle)*cylinder_area*(T(angle)-T_wall);
        H(angle) = mw(angle)*cp_total(angle)*T(angle);
        % finding the piston height and velocity
        piston_height(angle) = stroke - stroke*abs(sind(angle/2));
        deltay(angle) = stroke*sind(angle/2);
        pistonvel(angle) = deltay*alpha_per_second;
        % finding the work
        w(angle) = P_intake*cylinder_area*abs(deltay(angle));
    end
% calculating the new energy
    if angle == 1
        E(angle) = Eend - Q(angle) + W(angle) - H(angle);
    else
        E(angle) = E(angle-1) - Q(angle) + W(angle) - H(angle);
    end

end

mdottotal = sum(mdot);
Etotal = sum(E);
Qtotal = sum(Q);
Wtotal = sum(W);
Htotal = sum(H);
T_exhaust_pipe_total = sum(T_exhaust_pipe);
P_exhaust_pipe_total = sum(P_exhaust_pipe);
W_totaltotal = sum(W_total);
Q_totaltotal = sum(Q_total);