function [m,y,R_cw] = Counterweights(R, h, t, a, o, r, rho) % returns mass and center of mass radial-distance for counterweight
l = (2*R*sind(a/2) - 2*r)/2;
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

alpha = 0; %angle between counterweight and crankpin, zero for a flat twin
Rtop = h + abs(l)*tan(o) + R*(1-cosd(alpha/2));
Rside = sqrt((l+r)^2 + (h+abs(l)*tan(o))^2);
R_cw = max(Rtop,Rside); % furthest protrusion of counterweight from rotational axis

end