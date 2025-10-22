function [] = FiniteHeatMassLoss()
% Gas cycle heat release code with and w/o heat transfer
% data structure for engine parameters
clear();
thetas = -20; % start of heat release (deg)
thetad = 40; % duration of heat release (deg)
r =10;       % compression ratio
gamma = 1.4; % gas const
Q = 20.;    % dimensionless total heat release
h = 0.2;    % dimensionless ht coefficient
tw = 1.2;   % dimensionless cylinder wall temp
beta = 1.5; % dimensionless volume
a = 5;       % weibe parameter a 
n = 3;       % weibe exponent n
omega =209.4; % engine speed rad/s
c = 0.8;     % mass loss coeff

step=1;     % crankangle interval for calculation/plot
NN=360/step; % number of data points

theta = -180; % initial crankangle
thetae = theta + step; % final crankangle in step

% initialize results data structure
save.theta=zeros(NN,1);
save.vol=zeros(NN,1);   % volume 
save.press=zeros(NN,1); % pressure 
save.work=zeros(NN,1);  % work
save.heatloss=zeros(NN,1); % heat loss
save.mass=zeros(NN,1);  % mass left
fy=zeros(4,1); % vector for pressure, work, heat and mass loss
fy(1) = 1; % initial pressure (bar)
fy(4) = 1; % initial mass (-)

%for loop for pressure and work calculation
for i=1:NN,
[fy, vol] = integrate_ht(theta,thetae,fy);
% print values
% fprintf('%7.1f  %7.2f   %7.2f   %7.2f \n', theta,vol,fy(1),fy(2),fy(3));

% reset to next interval
theta = thetae;
thetae = theta+step;

save.theta(i)=theta; % put results in output vectors
save.vol(i)=vol;
save.press(i)=fy(1);
save.work(i)=fy(2);
save.heatloss(i)=fy(3);
save.mass(i)=fy(4);
end % end of pressure and work for loop

[pmax, id_max] = max(save.press(:,1)); % find max pressure
thmax=save.theta(id_max);              % and crank angle
ptdc=save.press(NN/2)/pmax;
w=save.work(NN,1);       % w is cumulative work
massloss =1- save.mass(NN,1);
eta=w/Q;                 % thermal efficiency
imep = eta*Q*(r/(r -1)); %imep/P1V1
eta_rat = eta/(1-r^(1-gamma));

% output overall results
fprintf(' Weibe Heat Release with Heat and Mass Loss  \n');
fprintf(' Theta_start =      %5.2f  \n', thetas);
fprintf(' Theta_dur =         %5.2f  \n', thetad);
fprintf(' P_max/P1 =          %5.2f  \n', pmax);
fprintf(' Theta @P_max =   %7.1f  \n',thmax);
fprintf(' P_tdc/P_max =       %5.2f  \n', ptdc);
fprintf(' Net Work/P1V1 =   %7.2f  \n', w);
fprintf(' Heat Loss/P1V1 =  %7.2f  \n', save.heatloss(NN,1));
fprintf(' Mass Loss/m =      %7.3f  \n',massloss );
fprintf(' Efficiency =        %5.3f   \n', eta);
fprintf(' Eff./Eff. Otto =    %5.3f  \n', eta_rat);
fprintf(' Imep/P1 =           %5.2f  \n', imep);

%plot results

plot(save.theta,save.press,'-','linewidth',2 )
set(gca, 'fontsize', 18,'linewidth',1.5);
xlabel('Crank Angle \theta (deg)','fontsize', 18)
ylabel('Pressure P (bar)','fontsize', 18)


figure();
plot(save.theta,save.work,'-',save.theta,save.heatloss,'--','linewidth',2 )
set(gca, 'Xlim',[-180 180],'fontsize', 18,'linewidth',1.5);
hleg1=legend('Work', 'Heat Loss','Location','NorthWest')
set(hleg1,'Box', 'off')
xlabel('Crank Angle \theta (deg)','fontsize', 18)
ylabel('Cumulative Work and Heat Loss','fontsize', 18)


function[fy,vol] = integrate_ht(theta,thetae,fy)
%  ode23 integration of the pressure differential equation
%  from theta to thetae with current values of fy as initial conditions

[tt, yy] = ode23(@rates, [theta thetae], fy);

% put last element of yy into fy vector
 for j=1:4
  fy(j) = yy(length(tt),j);
 end

 % pressure differential equation
    function [yprime] = rates(theta,fy) 
    vol=(1.+ (r -1)/2.*(1-cosd(theta)))/r;
    dvol=(r - 1)/2.*sind(theta)/r*pi/180.; %dvol/dtheta
    dx=0.;
        if(theta>thetas) % heat release >0
        dum1=(theta -thetas)/thetad;
        x=1-exp(-(a*dum1^n));
        dx=(1-x)*a*n*dum1^(n-1)/thetad; %dx/dthetha
        end
        
     term1= -gamma*fy(1)*dvol/vol;
     term3= h*(1. + beta*vol)*(fy(1)*vol/fy(4) - tw)*pi/180.;
     term2= (gamma-1)/vol*(Q*dx - term3);
     yprime(1,1)= term1 + term2 - gamma*c/omega*fy(1)*pi/180;
     yprime(2,1)= fy(1)*dvol;
     yprime(3,1)= term3;
     yprime(4,1)= -c*fy(4)/omega*pi/180;
    end %end of function rates
end % end of function integrate_ht
end % end of function HeatReleaseHeatTransfer
