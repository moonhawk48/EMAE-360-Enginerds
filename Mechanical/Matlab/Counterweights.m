function [m,y,R_top,R_side] = Counterweights(R, h, t, alph, del, r, rho) % returns mass and center of mass radial-distance for counterweight
l = R*sin(alph/2) - r;
lam = 0;
gam = 0;
if (h + abs(l)*tan(del)) <= r
    lam = 2*acos((h + abs(l)*tan(del))/r);
end
if h <= r
    gam = 2*acos(h/r);
end
A1 = (R^2)*(alph/2 - sin(alph/2)*cos(alph/2)) - (r^2)*(lam/2 - sin(lam/2)*cos(lam/2));
A2 = (2*r + l)*abs(l)*tan(del) - (r^2)*(gam/2 - sin(gam/2)*cos(gam/2)) + (r^2)*(lam/2 - sin(lam/2)*cos(lam/2));
A3 = 2*r*h - (r^2)*pi/2 + (r^2)*(gam/2 - sin(gam/2)*cos(gam/2));
A = A1 + A2 + A3;
m = A*t*rho/1000; % return mass in kg from rho given in g/cm^3


if (A1 == 0)
    y1 = 0;
else
    y1 = (1/A1)*((R^2)*(alph/2 - sin(alph/2)*cos(alph/2))*(h + abs(l)*tan(del) - R*cos(alph/2) + (((2*R*sin(alph/2))^3)/(12*(R^2)*(alph/2 - sin(alph/2)*cos(alph/2))))) - ((2/3)*(r^3)*(sin(lam/2))^3));
end
if (A2 == 0)
    y2 = 0;
else
    y2 = (1/A2)*((2*r + l)*abs(l)*tan(del)*(h + abs(l)*tan(del) - (abs(l)*tan(del)*(3*r + l))/(3*(2*r + l))) - ((2/3)*(r^3)*(sin(gam/2))^3 - (sin(lam/2))^3));
end
if (A3 == 0)
    y3 = 0;
else
    y3 = ((h^2)*r - ((r^3)*((2/3) - (2/3)*(sin(lam/2))^3)))/A3;
end

y = (y1*A1 + y2*A2 + y3*A3)/(A1+A2+A3);

R_top = h + abs(l)*tan(del) + R*(1-cos(alph/2));
R_side = sqrt((l+r)^2 + (h+abs(l)*tan(del))^2);

end