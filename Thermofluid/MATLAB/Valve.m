
function [mdot,choke] = Valve(stage,B,T_0,P_0,P_t,cp,cv)
%inlet valve function
%valve dimension array: See Heywood Ch6. 
%arbitrary valve geometry parameters

%inputs: open stage, bore size, inlet temp, inlet
%pressure, outlet pressure, cp, cv
libLoc = "C:\Users\Lucian\AppData\Roaming\CoolProp";

%Using Heywood terminology: P_0 = inlet pressure, P_t = outlet pressure

    %valve geometry variables
    % Valve Diameter, Seat angle, Seat Diameter,Port diameter, stem diameter
    Dv = 0.43*B;
    beta = deg2rad(45);
    Dm = 0.9*Dv;
    Dp = 1.1*Dv;
    Ds = 0.25*Dv;
    w = (Dv - Dm)/2;
    R = 8.341
    
    if stage == 1 %small opening
        Lv = w/(sin(beta)*cos(beta))
    elseif stage == 2 %partially open
            Lv = (sqrt((Dp^2-Ds^2)/(4*Dm)-w^2)+w*tan(beta))*(sin(beta)*cos(beta))/w
    elseif stage == 3 %fully open
            Lv = (sqrt((Dp^2-Ds^2)/(4*Dm)-w^2)+w*tan(beta))
    end
    Ac = pi*Dv*Lv;

    y = cp/cv;
    Cd = 0.8;   %discharge coefficient, experimentally determined
        %may need to vary based on inlet/outlet and choked vs unchoke
    

    %choked flow check
    if P_t/P_0 <= (2/(y+1))^(y/(y-1))
        choke = 1
        mdot = Cd*Ac*P_0/sqrt(R*T_0)*sqrt(y)*(2/(y+1))^(y+1)/(2*(y+1))
    else
        choke = 0
        mdot = (Cd*Ac*P_0/R*sqrt(T_0))*((P_t/P_0)^(1/y))*sqrt((2*y/(y-1))*(1-(P_t/P_0)^(y-(1/y))))
    end
    
end