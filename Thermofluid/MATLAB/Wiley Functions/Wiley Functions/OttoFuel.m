%program OttoFuel - computes const vol fuel air cycle

% first,isentropic compression from v1 to known v2 

%establish initial conditions at state 1
clear;
T1 = 350; %Kelvin
P1 = 101.3; %kPa
phi = 1.1; %equivalence ratio
f= .1; %residual fraction, 
rc=10; %compression ratio
    
fuel_id=2 ; %id: 1=Methane, 2=Gasoline, 3=Diesel, 4=Methanol, 5=Nitromethane
[ ~, ~, ~, ~, ~, ~, ~, ~, FS, ac ] = fuel( fuel_id, T1 );
% FS is stoichiometric fuel/air ratio
% ac is available energy of combustion kJ/kg

maxits = 50; 
tol =0.0001;

% call farg to get properties at 1
[y1, h1,u1, s1, v1, r, cp1, mw, dvdT, dvdP] = farg( T1, P1, phi, f, fuel_id );

%initial estimates of T2,P2
v2=v1/rc;
s2=s1;
cv1=cp1+ T1*(dvdT^2)/dvdP;
gam= cp1/cv1;
T2=T1*(v1/v2)^(gam-1.);
P2=P1*(v1/v2)^gam;

%do the iteration to find T2 and P2
for i2 = 1:maxits,
    [y2, h2,u2, s2, v2, r, cp2, mw, dvdT, dvdP] = farg( T2, P2, phi, f, fuel_id );
    f1=s1-s2;
    f2=v1/rc - v2;
    det= cp2*dvdP/T2 + dvdT^2;
    dt=(dvdP*f1 + dvdT*f2)/det;
    dp= (-dvdT*f1 + cp2/T2*f2)/det;
    %update T2 and P2
    T2=T2 + dt;
    P2=P2 + dp;
    %check for convergence
    %check for convergence
    if ( abs(dt)/T2 < tol && abs(dp)/P2 < tol )
        break;
    end
end

w12=-(u2-u1);% compression work

%combustion from 2-3 with v and u constant

%initial estimates of T3,P3 at state 3 
T3=3000;%Kelvin
P3=7000; % kPa

%do the iteration to find T3 and P3
for i3 = 1:maxits,
    [ierr,y3, h3, u3, s3, v3, r, cp3, mw, dvdT, dvdP] = ecp( T3, P3, phi, fuel_id );
    %fprintf(' \n combustion ierr= %6.2f  \n', ierr);
    f1= u2-u3;
    f2= v2-v3;
    det= cp3*dvdP+T3*dvdT^2;
    dt= (-f1*dvdP - f2*(dvdT+dvdP))/det;
    dp= ((P3*dvdT-cp3)*f2 + f1*dvdT)/det;
    %update T3 and P3  
    T3=T3 - dt;
    P3=P3 - dp;
    %check for convergence
    %check for convergence
    if ( abs(dt)/T3 < tol && abs(dp)/P3 < tol )
        break;
    end
end

% isentropic expansion from v3 to known v4 

%initial estimates of T4,P4
v4=v1;
cv3=cp3+ T3*(dvdT^2)/dvdP;
gam= cp3/cv3;
T4=T3*(v3/v4)^(gam-1.);
P4=P3*(v3/v4)^gam;
%T4=1500.;
%P4=500.;

%do the iteration to find T4 and P4
for i4 = 1:maxits,
    [ierr, y4, h4, u4, s4, v4, r, cp4, mw, dvdT, dvdP] = ecp( T4, P4, phi, fuel_id );
     %fprintf(' \n expansion ierr= %6.2f  \n', ierr);
    f1=s3-s4;
    f2=rc*v3 - v4;
    det= cp4*dvdP/T4 + dvdT^2;
    dt=(dvdP*f1 + dvdT*f2)/det;
    dp= (-dvdT*f1 + cp4/T4*f2)/det;
    %update T4 and P4
    T4=T4 + dt;
    P4=P4 + dp;
    %check for convergence
    if ( abs(dt)/T4 < tol && abs(dp)/P4 < tol )
        break;
    end
end
%compute cycle parameters
w=u1-u4;% net work
eta=w*(1+phi*FS)/phi/FS/(1.-f)/ac;
imep=w/(v1-v2); %imep

%output state and cycle parameters
fprintf(' \n Ottofuel input conditions: phi= %6.2f  fuel= %4d \n', phi, fuel_id);
fprintf(' State \t\t      1 \t    2 \t \t    3\t  \t     4 \n')
fprintf(' Pressure (kPa)=  %7.1f \t %7.1f \t %7.1f \t %7.1f \n',P1,P2,P3,P4);
fprintf(' Temperature (K)= %7.1f \t %7.1f \t %7.1f \t %7.1f  \n',T1,T2,T3,T4);
fprintf(' Enthalpy(kJ/kg)= %7.1f \t %7.1f \t %7.1f \t %7.1f  \n',h1,h2,h3,h4);
fprintf(' Int.Energy(kJ/kg)=%7.1f \t %7.1f \t %7.1f \t %7.1f  \n',u1,u2,u3,u4);
fprintf(' Volume (m^3/kg) =%7.3f \t %7.3f \t %7.3f \t %7.3f  \n', v1,v2,v3,v4);
fprintf(' Entropy(kJ/kg K)=%7.3f \t %7.3f \t %7.3f \t %7.3f  \n',s1,s2,s3,s4);
fprintf(' Cp (kJ/kg K) =   %7.3f \t %7.3f \t %7.3f \t %7.3f  \n \n', cp1,cp2,cp3,cp4);
fprintf(' Work (kJ/kg)= %7.1f \n', w);
fprintf(' Efficiency= %7.3f \n', eta);
fprintf(' Imep (kPa)= %7.1f \n \n', imep);
fprintf(' Iterations = \t \t %4d \t %4d \t %4d   \n', i2,i3,i4);

