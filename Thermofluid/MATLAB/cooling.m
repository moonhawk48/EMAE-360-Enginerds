
function [Tblock] = cylcooling(Twall,h1,L,cyldia,liner_thickness,block_thickness,Qin,ka,kb)


r1 = .5*cyldia;
r2 = r1+liner_thickness;
r3 = r2 + block_thickness;

rtot = (1/(2*pi*r1*L*h1)) + ((log(r2/r1))/(2*pi*ka*L)) + ((log(r3/r2))/(2*pi*kb*L));

theta = Qin*rtot;

Tblock = Twall-theta;


end