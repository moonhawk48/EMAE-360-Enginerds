function [Rmin, Kf, Km] = Counterweights(r,rcw,mrot,rho,Ycw,a, R,h,t,o)
%r is the crank radius
%rcw is the base counterweight radius around crank pin
%mrot is the rotational mass of crankshaft, accounting for webs, pins, and conrod bearing end
%rho is the density of counterweights
%Ycw = [1001] or [1111] depending on if there are counterweights in between or just outside

%a is distance between cranks midpoint-to-midpoint
w = (5000/60); %Rpm to Rps

theta = 0; %angle of first crank, I think should be 0
phi = 180; %phase between cranks based on firing order -- could also be 0
alpha = 0; %angle between counterweight and crankpin, I'm pretty sure this should be zero for a flat twin

CS = [sind(theta) sind(theta+phi)
      -a/2              a/2
      cosd(theta) cosd(theta+phi)];

CW = [-Ycw(1:2)*sind(theta-alpha), -Ycw(3:4)*sind(theta+phi+alpha)
      Ycw(1)*(CS(2,1)-a/4), Ycw(2)*(CS(2,1)+a/4), Ycw(3)*(CS(2,2)-a/4), Ycw(4)*(CS(2,2)+a/4)
      -Ycw(1:2)*cosd(theta-alpha), -Ycw(3:4)*cosd(theta+phi+alpha)];

l = (2*R*sind(a/2) - 2*r)/2;
[mcw,ycw] = mcw_def(R,h,t,a,o,rcw,rho,l);

Fcw = sqrt((sum(mcw*(w^2)*ycw*CW(1,1:4)))^2 + (sum(mcw*(w^2)*ycw*CW(3,1:4)))^2);
Fcs = sqrt((sum(mrot*(w^2)*r*CS(1,1:2)))^2 + (sum(mrot*(w^2)*r*CS(1,1:2)))^2);
if (Fcs == 0 && Fcw == 0)
    Kf = 100; %if both forces are internally balanced, default to 100% balanced
else
    Kf = 100*Fcw/Fcs; % centrifugal(?) forces balanced by counterweights (should be 0/0 for boxer since inherently balanced)
end
Mcw = sqrt((sum(mcw*(w^2)*ycw*CW(1,1:4).*CW(2,1:4)))^2 + (sum(mcw*(w^2)*ycw*CW(3,1:4).*CW(2,1:4)))^2);
Mcs = sqrt((sum(mrot*(w^2)*r*CS(1,1:2).*CS(2,1:2)))^2 + (sum(mrot*(w^2)*r*CS(3,1:2).*CS(2,1:2)))^2);
Km = 100*Mcw/Mcs;  % moment forces balanced by counterweights

Rtop = h + abs(l)*tan(o) + R*(1-cosd(alpha/2));
Rside = sqrt((l+rcw)^2 + (h+abs(l)*tan(o))^2);
Rmin = max(Rtop,Rside);
end

function [m,y] = mcw_def(R, h, t, a, o, r,rho,l) % returns mass and center of mass radial-distance for counterweight

lam = 0;
gam = 0;
if (h + abs(l)*tan(o)) <= r
    lam = 2*acosd((h + abs(l)*tan(o))/r);
end
if h <= r
    gam = 2*acosd(h/r);
end
A1 = (R^2)*(a/2 - sind(a/2)*cosd(a/2)) - (r^2)*(lam/2 - sind(lam/2)*cosd(lam/2));
A2 = (2*r + l)*abs(l)*tan(o) - (r^2)*(gam/2 - sind(gam/2)*cosd(gam/2)) + (r^2)*(lam/2 - sind(lam/2)*cosd(lam/2));
A3 = 2*r*h - (r^2)*pi/2 + (r^2)*(gam/2 - sind(gam/2)*cosd(gam/2));
A = A1 + A2 + A3;
m = A*t*rho;


y1 = (1/A1)*((R^2)*(a/2 - sind(a/2)*cosd(a/2))*(h + abs(l)*tan(o) - R*cosd(a/2) + ((2*R*sind(a/2)^3)/(12*(R^2)*(a/2 - sind(a/2)*cosd(a/2))))) - ((2/3)*(r^3)*(sind(lam/2))^3));
y2 = (1/A2)*((2*r + l)*abs(l)*tan(o)*(h + abs(l)*tan(o) - (abs(l)*tan(o)*(3*r + l))/(3*(2*r + l))));
y3 = ((h^2)*r - ((r^3)*((2/3) - (2/3)*(sind(lam/2))^3)))/A3;

y = (y1*A1 + y2*A2 + y3*A3)/(A1+A2+A3);

end