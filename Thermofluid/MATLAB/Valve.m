%inlet valve function
%valve dimension array: See Heywood Ch6. 
%arbitrary valve geometry parameters
% Valve Diameter, Seat angle, Seat Diameter,Port diameter, stem diameter
%dim = [0.43*B,45,0.9*Dv,Dv*1.1,.25*Dv]
%inputs: open stage, bore size, valve geomewtry array, inlet temp, inlet
%pressure, outlet pressure, fluid type, Coolprop Library location

function [mdot,choke] = Valve(open,B,T_0,P_0,P_t)
%Using Heywood terminology: P_0 = inlet pressure, P_t = outlet pressure

    %valve geometry variables
    Dv = 0.43*B;
    beta = deg2rad(45);
    Dm = 0.9*Dv;
    Dp = 1.1*Dv;
    Ds = 0.25*Dv;
    w = (Dv - Dm)/2;
    
    if open == 1 %partially open
        Lv = w/(sin(beta)*cos(beta));
    elseif open == 2 %fully open
        Lv = sqrt(((Dp^2-Ds^2)/(4*Dm))^2-w^2)+w*tan(beta)
    end

    Ac = pi*Dv*Lv;

    % cp = getFluidProperty(libLoc,'CPMASS','T',T_0,'P',P_0,fluid);
    % cv = getFluidProperty(libLoc,'CVMASS','T',T_0,'P',P_0,fluid);
    % y = cp/cv;
    
    R = 8.341;
    %Egnell function
    y = 1.38-0.2*exp(-900/T_0);
    Cd = 0.8;   %discharge coefficient, experimentally determined
    %may need to vary based on inlet/outlet and choked vs unchoked

    %choked flow check
    if P_t/P_0 <= (2/(y+1))^(y/(y-1))
        choke = 1;
        mdot = Cd*Ac*P_0/sqrt(R*T_0)*sqrt(y)*(2/(y+1))^(y+1)/(2*(y+1));
    else
        choke = 0;
        mdot = (Cd*Ac*P_0/(R*sqrt(T_0)))*((P_t/P_0)^(1/y))*sqrt((2*y/(y-1))*(1-(P_t/P_0)^(y-(1/y))));
    end
    
end