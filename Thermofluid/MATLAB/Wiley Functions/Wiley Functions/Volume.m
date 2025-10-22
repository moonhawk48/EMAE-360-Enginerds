function []=Volume()
% this program computes and plots the exact 
% and approximate cylinder volume
clear();
r = 8; % compression ratio
s = 100;   % stroke (mm)
len= 150; %connecting rod length (mm)
ep=s/(2*len);
theta=-180:1:180; %crankangle theta vector
ys1=(1-cosd(theta))/2; %approx y/s
ys2= ys1+ (1-(1- ep^2*sind(theta).^2).^(1/2))/(2*ep); %exact y/s
vol1 = 1+(r-1)*ys1; %approx volume
vol2= 1+(r-1)*ys2;  % exact volume
diff = abs(vol1-vol2)./vol1 * 100;
[diffmax, id_max] = max(diff);
thmax=theta(id_max);

%plot results
fprintf('Max Diff (percent) = %5.2f  at theta (deg) = %5.2f \n', diffmax,thmax);
figure;
plot(theta,vol1,'--',theta,vol2,'-','linewidth',2);
set(gca,'Xlim',[-180 180],'Ylim',[0 r],'fontsize',18,'linewidth',2);
xlabel('Crank Angle (deg)','fontsize', 18);
ylabel('Dim. Cylinder Volume','fontsize', 18);
legend('Approx. Volume', 'Exact Volume','Location', 'North');
figure;
plot(theta,diff,'linewidth',2);
set(gca,'Xlim',[-180 180],'fontsize',18,'linewidth',2);
xlabel('Crank Angle (deg)','fontsize', 18);
ylabel('Dim. Error (%)','fontsize', 18);

end