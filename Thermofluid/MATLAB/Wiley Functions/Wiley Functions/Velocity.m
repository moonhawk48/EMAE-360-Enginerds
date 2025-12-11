function []=Velocity()
% this program computes and plots the 
% dimensional and dimensionless piston velocity
clear();
N=2000; %rev/min
s = 0.10;   % stroke (m)
len= 0.150; %connecting rod length (m)
ep=s/(2*len);

theta=0:1:180; %crankangle theta vector
term1=pi/2*sind(theta);
term2= (1+(ep*cosd(theta))./(1 - ep^2*sind(theta).^2).^(1/2)); %exact y/s
Upbar=term1.*term2; %dimensionless velocity
Up=Upbar*2*N*s/60; % m/s
[umax, id_max] = max(Up);% maximum velocity
thmax=theta(id_max); % crankangle for max velocity

%plot results
fprintf( 'U max (m/s) = %5.2f  at theta (deg) = %5.2f \n', umax,thmax);
figure;
plot(theta,Up,'linewidth',2);
set(gca,'Xlim',[0 180],'fontsize',18,'linewidth',2);
xlabel('Crank Angle (deg)','fontsize', 18);
ylabel(' Piston Velocity (m/s)','fontsize', 18);
figure;
plot(theta,Upbar,'linewidth',2);
set(gca,'Xlim',[0 180],'Ylim',[0 1.8],'fontsize',18,'linewidth',2);
xlabel('Crank Angle (deg)','fontsize', 18);
ylabel('Dimensionless Piston Velocity','fontsize', 18);
end