% Four stroke Otto cycle model
% input parameters
Ti = 300; % inlet temperature, K
Pi = 50; % inlet pressure, kPa
Pe = 100; % exhaust pressure, kPa
r = 10; % compression ratio
qin = 2500; % heat input, kJ/kg (mixture)
gamma = 1.3; % ideal gas specific heat ratio
R = 0.287; % gas constant kJ/kg K
cv= R/(gamma-1); %const vol specific heat, kJ/kg K

f=0.05;% guess value of residual fraction f
Te = 1000; % guess value of exhaust temp, K
tol=0.0001; % tolerance for convergence
err = 2*tol; %error initialization
gam=(gamma -1)/gamma;

while (err > tol) %while loop for cycle calc
%intake stroke
T1=(1-f)*Ti + f*(1 - (1- Pi/Pe)*gam)*Te;
P1=Pi;
%isentropic compression
P2=P1*r^gamma;
T2=T1*r^(gamma-1);
%const v heat addition
T3=T2 + qin*(1-f)/cv;
P3=P2*(T3/T2);
%isentropic expansion
P4=P3*(1/r)^gamma;
T4=T3*(1/r)^(gamma-1);
%isentropic blowdown
P5=Pe;
T5=T4*(P4/Pe)^(-gam);
%const p exhaust stroke
Te=T5;
fnew=(1/r)*(Pe/P4)^(1/gamma); %new residual fraction
err=abs(fnew-f)/fnew;
f=fnew;
end

%cycle parameters
eta= 1 - 1/r^(gamma-1);
imep = P1*(qin*(1-f)/(R*T1))/(1-(1/r))*eta;
pmep=Pe-Pi;
etanet= eta*(1-pmep/imep);
imepnet= (imep-pmep)/100.;
voleff=1-(Pe/Pi -1)/(gamma*(r-1));

%output calcs
fprintf(' \nFour Stroke Otto Cycle  \n')
fprintf('State              1           2           3         4 \n');
fprintf('Pressure (kPa):   %6.1f     %6.1f     %6.1f    %6.1f \n',P1,P2,P3,P4);
fprintf('Temperature (K):  %6.1f     %6.1f     %6.1f    %6.1f \n',T1,T2,T3,T4);
fprintf('Ideal Thermal Eff.= %6.3f         Net Thermal Eff.=     %6.3f \n',eta, etanet);
fprintf('Exhaust Temp. (K)=  %6.1f         Volumetric Eff.=      %6.2f \n',Te, voleff);
fprintf('Residual Fraction   %6.3f         Net Imep (bar)=      %6.2f \n',f, imepnet);


