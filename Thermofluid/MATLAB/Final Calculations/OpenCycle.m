function [out]=OpenCycle(speed)
%calculates the Open cycle portion of the Otto cycle

clear
% engine speed (rpm)
%speed= [800,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000];
workRPM = [];
speed = [800,5000,6000];
% relative valve timing angles (deg)
io_btc=15; % intake open before top center
ic_abc=0; % intake closed after bottom center
eo_bbc=30; % exhaust open before bottom center
ec_atc=15; % exhaust closed after top center
%cylinder specs.
s=0.0864; % stroke (m)
a=s/2; % crank throw radius (m)
b=0.108; % bore (m)
R=3.*a; % rod length (m)
cr=10; % compression ratio
%intake & exhaust valve specifications
di = 0.05184;    % diameter of intake seat (m)
Li = 0.01803;    % maximum inlet valve lift (m)
de = 0.04644;    % diameter of exhaust seat (m)
Le = 0.01615;    % maximum exhaust valve lift (m)
Pe=105; % exhaust pressure (kPa)
Pi=100; % intake pressure (kPa)
vd=pi/4*(b^2)*s; % displacement volume per cylinder (m^3)
vbdc= vd/(1-1/cr); % volume at bottom dead center
vc=vd/(cr-1); % clearance volume (m^3)
%Energy release specs.
Qbar=20; % dimensionless heat release Qin/P1V1
Qin=Qbar*Pi*vbdc; % total heat release(kJ/cycle)
a1=5; % wiebe form factor
n=3; % wiebe efficiency factor
thetas=-35; % start of heat release (deg)
thetad=60; % heat release duration (deg)
Rs=0.287; % air gas constant (kJ/kg-K)
Ri=Rs; Re=Rs; % intake and exhaust gas constant(KJ/kg-K)
ka=1.4; % specific heat ratio
ki=ka; ke=ka; % specific heat ratio for intake and exhaust gas
%convert to absolute crankangle degrees
IVO=360-io_btc;
IVC=540+ic_abc;
EVO=180-eo_bbc;
EVC=360+ec_atc;
thetadi=IVC-IVO; % duration of valve opening (deg)
thetade=EVC-EVO; % duration of valve opening (deg)
Lthetai=0; % initial condition for intake valve lift
Lthetae=0; % initial condition for exhaust valve lift
start=IVC-720; % start of simulation
dtheta=1; % theta change per step
NN=360/dtheta; % number of crank angle data points in calc.
NIVO= (IVC-IVO)/dtheta; % number of crank angle data points with intake valve open
NEVO= (EVC-EVO)/dtheta; % number of crank angle data points with exhaust valve open
%initalize and allocate vector space
wtheta=zeros(NN,1); %crank angle vector (deg)
wPc=zeros(NN,2); %pressure (kPa)
wvol=zeros(NN,1); %volume (m^3)
wtem=zeros(NN,2); %cylinder temperature (K)
wMi=zeros(NN,2); %mass in cylinder (kg)
wMii=zeros(NN,2); %inlet mass in/out of cylinder (kg)
wwork=zeros(NN,2); %differential work (kJ)
wLthetai=zeros(NN,1); %intake valve lift (mm)
wLthetae=zeros(NN,1); %exhaust vave lift (mm)
iPc=zeros(NIVO,2); %intake valve open pressure (kPa)
itheta=zeros(NIVO,1); %intake valve open crank angle vector (deg)
idmi=zeros(NIVO,2); %intake valve open inlet mass vector (kg/deg)
ePc=zeros(NEVO,2); %exhaust valve open pressure (kPa)
etheta=zeros(NEVO,1); %exhaust valve open crank angle vector (deg)
edme=zeros(NEVO,2); %exhaust valve open inlet mass vector (kg/deg)
% main for loop for the two engine speeds
for k=1:length(speed)
    N=speed(k); %engine speed (rpm)
    omega=360*N/60; %angular velocity (deg/s)
    Up=4*N*a/60; %average piston speed (m/s)
    %initial conditions
    Mi=0.000866; %initial guess of mass (kg) in cylinder
    Pc=100; %initial guess of cylinder pressure Pc (kPa)
    Tc=344; %initial guess of cylinder temperature Tc (K)
    Ti=320; %initial intake temp.(K)
    Te=600; %inital exhaust temp.(K)
    Tw=500; %cylinder wall temp. (K)
    C1=Ti;
    C2=Te;
    %the iteration loop to find converged steady state
    for iteration=1:8 % usually <5 iterations
        i=1; %reset step counter to 1
        j=1; %counter of intake
        jj=1; %and exhaust steps
        Mii=0; %initial mass into cylinder through intake valve
        Tcold=Tc; %save Tc(ivc) from previous iteration
        Pcold=Pc; %save Pc(ivc) from previous iteration
        Miold=Mi; %save Mi(ivc) from previous iteration
        tol=0.0001;
        % guess cylinder conditions at start=ivc
        [Vivc,dvol,yc]=geofun (R, b, a, start, vc, vd);
        Mtr=Pc*Vivc/(Rs*Tc); %trapped mass at ivc
        Mi=Mtr; %mass in cylinder at ivc is trapped mass
        fy(1)=Pc;
        fy(2)=Mi; % guess values of P, mass, and work at ivc
        fy(3)=0;
        % For loop to evaluate heat release from ivc to evo
        for theta=start:EVO-1 %
            thetae=theta+dtheta;
            % integrate closed energy eqn for pressure fy(1),and cylinder work fy(3)
            [fy,vol,Tc]=integrate2(theta, thetae, fy);
            %copy results to output vectors
            wtheta(i)=theta; %whole 4-stroke theta (deg)
            wPc(i,k)=fy(1); %pressure (kPa)
            wvol(i)=vol; %volume (m^3)
            wtem(i,k)=Tc; %cylinder temperature (K)
            wMi(i,k)=Mi; %mass in cylinder (kg)
            wwork(i)=fy(3); %differential work (kJ) %advance to next crankangle
            i=i+1;
        end %end of for loop to evaluate heat release from ivc to evo
        % For loop from EVO to IVC to compute intake and exhaust strokes
        fy(1)=Pc;
        fy(2)=Mi; % initial values at EVO
        for theta=EVO:IVC
            thetae=theta+dtheta;
            % integrate open energy eqn for pressure fy(1), and cylinder mass fy(2)
            [fy,vol,dmi,dme,Tc]=integrate(theta, thetae, fy);
            % copy results to output vectors
            wtheta(i)=theta; %whole 4-stroke theta (deg)
            wPc(i,k)=fy(1); %pressure (kPa)
            wvol(i)=vol; %volume (m^3)
            wtem(i,k)=Tc; %cylinder temperature (K)
            wMi(i,k)=fy(2); %mass in cylinder (kg)
            wMii(i,k)=Mii; %inlet mass into cylinder (kg)
            wwork(i,k)=fy(3); % differential work (kJ)
            wLthetai(i)=Lthetai*1000; %intake valve lift (mm)
            wLthetae(i)=Lthetae*1000; %exhaust vave lift (mm)
            Mii=Mii+dmi*dtheta; % add incremental mass into cylinder through intake valve (kg)
            % for volumetric efficiency calc.
            % advance to next crank angle interval
            i=i+1;
            if theta==EVC % save i index at EVC
                ievc=i;
            end
            if theta>=IVO %save data for open intake valve
                iPc(j,k)=Pc;
                itheta(j)=theta;
                idmi(j,k)=dmi;
                j=j+1;
            end %end of if theta>=IVO
            if theta>=EVO && theta<=EVC %save data for open exhaust valve
                ePc(jj,k)=Pc;
                etheta(jj)=theta;
                edme(jj,k)=dme;
                jj=jj+1;
            end %end of if theta>=IVO
            if theta==IVO %save data at ivo - intake valve open
                Pcivo=Pc;
                Vivo=vol;
            end %end of if theta==IVO
        end %end of for loop of intake and exhaust process
        % calculate overall cycle parameters
        rho=Pi/(Ri*Ti); %reference density for volumetric efficiency
        Ev=Mii/(rho*vd); %volumetric efficiency
        f=wMi(ievc,k)/wMi(1,k); %residual fraction
        w=wwork(IVC,k); %cumulative work from ivc to ivc
        eta = w/Qin; %thermal efficiency
        imep= w/vd; % imep
        % test for convergence
        if abs(Tcold-Tc)/Tcold<=tol && abs(Pcold-Pc)/Pcold<=tol && abs(Miold-Mi)/Miold<=tol
            break;
        end %end of if condition to break out of iteration for loop
    end %end of the iteration loop for steady state convergence
    fprintf('Speed (rpm) = %6.0f iteration # %8.0f IVC Temperature (K) = %6.2f \n',N, iteration, Tc);
    fprintf('Trapped Mass (kg)= %10.6f iterated IVC Pressure (kPa) = %6.2f \n',Mtr, Pc);
    fprintf('Volumetric efficiency= %9.3f Residual fraction = %6.3f \n',Ev, f);
    fprintf('Thermal efficiency= %6.3f IMEP (kPa) = %9.2f \n',eta, imep);
    disp(max(wtem(:,k)))
    disp(max(wPc(:,k)))
end

% end of main engine speed k loop
% set up horizontal lines for plots
xe1= [90 450];
xe2= [330 600];
ye= [0 0];
xe3=[300 600];
yp=[Pi Pi];
% indices for valve plot
evoi=EVO-start+1;
evci=EVC-start;
ivoi=IVO-start+1;
ivci=IVC-start;
%outputs
out = [wtheta,wPc(:,1),wtem(:,1)];

%plotting

% %valve lift plot
% figure;plot(wtheta(ivoi:ivci), wLthetai(ivoi:ivci),'b','linewidth',2);
% hold on % put another line on plot
% plot(wtheta(evoi:evci), wLthetae(evoi:evci),'r','linewidth',2);
% plot([180 180], get(gca,'ylim'),'k','linewidth',1);
% plot([360 360], get(gca,'ylim'),'k','linewidth',1);
% plot([540 540], get(gca,'ylim'),'k','linewidth',1);
% hold off
% xlabel('Crank Angle (deg)','fontsize', 22);
% ylabel('Valve Lift (mm)','fontsize', 22);
% title('Valve Lift vs. CA','fontsize', 22);
% legend('Intake', 'Exhaust','location','northwest');
% grid on
% set(gca, 'YLim',[0 15])
% set(gca, 'fontsize', 20,'linewidth',1.5);
% set(gca, 'XTick',[ 90 180 270 360 450 540],'XMinorTick','on')
% % Overall pressure plot
% figure;
% plot(wtheta, wPc(:,1),'b','linewidth',2);
% hold on
% plot(wtheta, wPc(:,2),'r','linewidth',2);
% P_ex = [wtheta(evoi:evci),wPc((evoi:evci),:)];
% writematrix(P_ex,'Exhaust_pressures')
%
% plot([0 0], get(gca,'ylim'),'k','linewidth',1);
% hold off
% legend([num2str(speed(1)) 'rpm'],[ num2str(speed(2)) 'rpm'],'location', 'northwest');
% grid on
% xlabel('Crank Angle (deg)','fontsize', 22);
% ylabel('Cylinder Pressure (Kpa)','fontsize', 22);
% title('Cylinder Pressure vs. CA','fontsize', 22);
% set(gca, 'fontsize', 20,'linewidth',1.5);
% set(gca, 'XTick',[-180 -90 0 90 180, 270, 360, 540],'XMinorTick','on')
% % cylinder temperature plot
% figure;
% plot(wtheta(), wtem(:,1),'b','linewidth',2);
% hold on
% plot(wtheta(), wtem(:,2),'r','linewidth',2);
% plot([0 0], get(gca,'ylim'),'k','linewidth',1);
% hold off
% legend([num2str(speed(1)) 'rpm'],[ num2str(speed(2)) 'rpm'],'location', 'northwest');
% grid on
% xlabel('Crank Angle (deg)','fontsize', 22);
% ylabel('Cylinder Temperature (K)','fontsize', 22);
% title('Cylinder Temperature vs. CA','fontsize', 22);
% set(gca, 'fontsize', 20,'linewidth',1.5);
% set(gca, 'XTick',[-180 -90 0 90 180, 270, 360, 540],'XMinorTick','on')
% % overall mass plot
% figure;
% plot(wtheta, wMi(:,1),'b','linewidth',2);
% hold on
% plot(wtheta, wMi(:,2),'r','linewidth',2);
% plot([180 180], get(gca,'ylim'),'k','linewidth',1);
% plot([360 360], get(gca,'ylim'),'k','linewidth',1);
% plot([540 540], get(gca,'ylim'),'k','linewidth',1);
% hold off
% xlabel('Crank Angle (deg)','fontsize', 22);
% ylabel('Net Cylinder Mass (kg)','fontsize', 22);
% title('Net Cylinder Mass','fontsize', 22);
% legend([num2str(speed(1)) ' rpm'],[num2str(speed(2)) ' rpm'],'location', 'southwest');
% grid on
% set(gca, 'fontsize', 20,'linewidth',1.5);
% set(gca, 'XTick',[-180 0 180 360 540],'XMinorTick','on')
% % inlet mass flow plot
% figure;
% plot(itheta, idmi(:,1),'b','linewidth',2);
% hold on
% plot(itheta, idmi(:,2),'r','linewidth',2);
% plot(xe2, ye,'k','linewidth',1);
% plot([360 360], get(gca,'ylim'),'k','linewidth',1);
% plot([540 540], get(gca,'ylim'),'k','linewidth',1);
% hold off
% legend([num2str(speed(1)) ' rpm'],[ num2str(speed(2)) ' rpm'], 'location','northeast');
% grid on
% xlabel('Crank Angle (deg)','fontsize', 22);
% ylabel('Inlet Mass Flowrate (kg/deg)','fontsize', 22);
% title('Inlet Mass Flowrate vs. CA','fontsize', 22);
% set(gca, 'fontsize', 20,'linewidth',1.5);
% set(gca, 'XLim',[330 600])
% set(gca, 'XTick',[360 450 540 630],'XMinorTick','on')
% %set(gca, 'XTick',[-180 0 180 360 540],'XMinorTick','on')
% % exhaust mass flow plot
% figure;
% plot(etheta, edme(:,1),'b','linewidth',2);
% hold on
% plot(etheta, edme(:,2),'r','linewidth',2);
% plot(xe1, ye,'k','linewidth',1);
% plot([360 360], get(gca,'ylim'),'k','linewidth',1);
% plot([180 180], get(gca,'ylim'),'k','linewidth',1);
% hold off
% legend([num2str(speed(1)) ' rpm'],[ num2str(speed(2)) ' rpm'], 'location','southeast');
% grid on
% xlabel('Crank Angle (deg)','fontsize', 22);
% ylabel('Exhaust Mass Flowrate (kg/deg)','fontsize', 22);
% title('Exhaust Mass Flowrate vs. CA','fontsize', 22);
% set(gca, 'XLim',[90 450])
% set(gca, 'fontsize', 20,'linewidth',1.5);
% set(gca, 'XTick',[90 180 270 360 450],'XMinorTick','on')
% % intake stroke cylinder pressure plot
% figure;
% plot(itheta, iPc(:,1),'b','linewidth',2);
% hold on
% plot(itheta, iPc(:,2),'r','linewidth',2);
% plot(xe3, yp,'k','linewidth',1);
% plot([360 360], get(gca,'ylim'),'k','linewidth',1);
% plot([540 540], get(gca,'ylim'),'k','linewidth',1);
% hold off
% legend([num2str(speed(1)) ' rpm'],[ num2str(speed(2)) ' rpm'], 'location','northwest');
% grid on
% xlabel('Crank Angle (deg)','fontsize', 22);
% ylabel('Cylinder Pressure (kPa)','fontsize', 22);
% title('Cylinder Pressure During Intake Process vs. CA','fontsize', 22);
% set(gca, 'fontsize', 20,'linewidth',1.5);
% set(gca, 'XLim',[300 600])
% set(gca, 'XTick',[360 450 540 630],'XMinorTick','on')
% %set(gca, 'XTickLabel',{'270'; '360'; '450';'540';'630'})


%engine geometric function to evaluate instantaneous engine volume
    function [vol,dvol,yc]=geofun (R, b, a, theta, vc, vd)
        y=a*cosd(theta)+sqrt(R^2-(a^2)*(sind(theta))^2); %evaluate piston position
        vol=vc+(R+a-y)*pi*b^2/4; %evaluate cylinder volume (m^3)
        dvol=vd/2*sind(theta)*(1+cosd(theta)/ sqrt(((R/a)^2-(sind(theta))^2)))*pi/180; %dv/dtheta
        yc=(R+a-y)+vc/(pi/4*b^2); %instantaneous piston location
    end
% function to integrate closed‐system energy eqn during combustion
    function[fy,vol,Tc]=integrate2(theta, thetae, fy)
        % integrate from theta to thetae with pointer to function rates2
        % and current values of vector fy as initial conditions
        [tt,yy] = ode23(@rates2, [theta thetae], fy);
        % put last elements of yy into fy vector
        for kk=1:3
            fy(kk)=yy(length(tt),kk);
        end
        % energy equation for dP/dtheta
        function [yprime] = rates2(theta, fy)
            y=a*cosd(theta)+sqrt(R^2-(a^2)*(sind(theta))^2); %evaluate piston position
            vol=vc+(R+a-y)*pi*b^2/4; %evaluate cylinder volume (m^3)
            dvol=vd/2*sind(theta)*(1+cosd(theta)/sqrt(((R/a)^2- (sind(theta))^2)))*pi/180; %dv/dtheta
            Pc=fy(1); %current cylinder pressure
            Tc=Pc*vol/(Mi*Rs); %cylinder temperature from ideal gas equation
            %compute heat loss dQw with Woschni correlation
            %T0=293; %inlet temp. (K)
            T0=Tcold; %reference gas temp at IVC
            Tg=Pc*vol/(Rs*Mtr); %instantaneous gas temp. (K)
            Pm=Pi*(Vivc^ka)/vol^ka; %motored pressure (Kpa) P(ivc)~ Pi
            deltaPc=Pc-Pm; %instantaneous pressure rise (kPa)
            u=2.28*Up+0.00324*T0*vd*deltaPc/(Vivc*Pi); %instantaneous gas speed (m/s)
            hg=3.26*(Pc^0.8)*(u^0.8)*(b^-0.2)*(Tg^-0.55); %heat transfer coefficient (W/m^2-K)
            Aw=pi*b*yc + pi/2*(b^2); % heat transfer area (m^2)
            dQw=hg*Aw*(Tg-Tw)/(360*N/60)/1000; %dQw/dtheta
            %compute combustion energy release
            dXb=0; Xb=0; % set heat release fractions to zero
            if theta>thetas % if theta>thetas, evaluate dXb by wiebe function
                Xb=1-exp(-a1*((theta-thetas)/thetad)^n); %heat release fraction Xb
                dXb=n*a1*(1-Xb)*(((theta-thetas)/thetad)^(n-1))/thetad; % dXb/dtheta
            end
            dQin=Qin*dXb; %dQin/dtheta
            term1=-ka*Pc*dvol/vol;
            term2= (ka-1)*(dQin-dQw)/vol;
            yprime(1,1)=term1+term2; %dP/dtheta
            yprime(2,1)= Mi; %mass in cylinder
            yprime(3,1)=Pc*dvol; % d work/dtheta
        end %end of function rates2
    end % end of function integrate2
% function to integrate open‐system energy eqn during valve events
    function [fy, vol,dmi, dme, Tc] = integrate(theta,thetae,fy)
        % integrate from theta to thetae with pointer to function rates
        % and current values of fy as initial conditions
        [tt,yy] = ode23(@rates, [theta thetae], fy);
        % put last elements of yy into fy vector
        for kk=1:3
            fy(kk)=yy(length(tt),kk);
        end
        % energy equation for dP/dtheta
        function [yprime] = rates(theta, fy)
            y=a*cosd(theta)+sqrt(R^2-(a^2)*(sind(theta))^2); %evaluate piston position
            vol=vc+(R+a-y)*pi*b^2/4; %evaluate instaneous cylinder volume (m^3)
            dvol=vd/2*sind(theta)*(1+cosd(theta)/sqrt(((R/a)^2- (sind(theta))^2)))*pi/180; %dv/dtheta
            Pc=fy(1);
            Mi=fy(2); %current values of pressure and cylinder mass
            Tc=Pc*vol/(Mi*Rs); %cylinder temperature
            %exhaust conditions
            Ce=sqrt(ke*Rs*1000*Te); %speed of sound in exhaust pipe condition (m/s)
            Cc=sqrt(ke*Rs*1000*Tc); %speed of sound in cylinder condition (m/s)
            rhoe=Pe/(Rs*Te); %density in exhaust pipe (kg/m^3)
            rhoc=Pc/(Rs*Tc); %density in cyinder (kg/m^3)
            Cre=((ke+1)/2)^(ke/(ke-1)); %critical pressure ratio
            Lthetae=Le*(1- cosd(360*(theta-EVO)/thetade))/2; % exhaust valve lift
            Cde=0.650+1.059*(Lthetae/de)-3.70*(Lthetae/de)^2; %discharge coefficient
            Ae=Cde*pi*de*Lthetae; %effective exhaust valve area
            % intake conditions
            Ci=sqrt(ki*Rs*1000*Ti); %speed of sound in intake pipe condition (m/s)
            Cc=sqrt(ki*Rs*1000*Tc); %speed of sound in cylinder condition (m/s)
            rhoi=Pi/(Rs*Ti); %density in intake pipe (Kg/m^3)
            Cri=((ki+1)/2)^(ki/(ki-1)); %critical pressure ratio
            Lthetai=Li*(1-cosd(360*(theta-IVO)/thetadi))/2; %valve lift (m)
            Cdi=0.850-0.860*(Lthetai/di)-0.798*(Lthetai/di)^2; %discharge coefficient
            Ai=Cdi*pi*di*Lthetai; %effective intake valve area (m^2)
            % get mass flow past valves
            % only exhaust valve open
            if theta>=EVO && theta<IVO
                if Pc/Pe>=Cre %outward choked flow
                    Fe=-(2/(ke+1))^((ke+1)/2/(ke-1));
                    dme=Ae*Cc*rhoc*Fe/omega;
                    C1=0;
                    C2=Tc;
                elseif Pc/Pe<Cre && Pc/Pe>=1 %outward nonchoked flow
                    Fe=-sqrt(2/(ke-1)*((Pe/Pc)^(2/ke)-(Pe/Pc)^((ke+1)/ke)));
                    dme=Ae*Cc*rhoc*Fe/omega;
                    C1=0;
                    C2=Tc;
                elseif Pe/Pc>=Cre %inward choked flow
                    Fe=(2/(ke+1))^((ke+1)/2/(ke-1));
                    dme=Ae*Ce*rhoe*Fe/omega; %exhaust valve mass flow rate dm/dtheta (kg/deg)
                    C1=0; %two parameters in cylinder pressure change equation
                    C2=Te;
                else %inward nonchoked flow
                    Fe=sqrt(2/(ke-1)*((Pc/Pe)^(2/ke)-(Pc/Pe)^((ke+1)/ke)));
                    dme=Ae*Ce*rhoe*Fe/omega;
                    C1=0;
                    C2=Te;
                end % end Pe/Pc check
                dmi=0; % no mass flow rate through intake valve
            end % end of exhaust flow calc
            % overlap with both valves open
            if theta>=IVO && theta<=EVC
                if Pc/Pe>=Cre %outward exhaust choked flow
                    Fe=-(2/(ke+1))^((ke+1)/2/(ke-1));
                    dme=Ae*Cc*rhoc*Fe/omega;
                    C2=Tc;
                elseif Pc/Pe<Cre && Pc/Pe>=1 %outward nonchoked flow
                    Fe=-sqrt(2/(ke-1)*((Pe/Pc)^(2/ke)-(Pe/Pc)^((ke+1)/ke)));
                    dme=Ae*Cc*rhoc*Fe/omega;
                    C2=Tc;
                elseif Pe/Pc>=Cre %choked inward flow
                    Fe=(2/(ke+1))^((ke+1)/2/(ke-1));
                    dme=Ae*Ce*rhoe*Fe/omega; %exhaust valve mass flow rate dm/dtheta (kg/deg)
                    C2=Te;
                else %inward nonchoked flow
                    Fe=sqrt(2/(ke-1)*((Pc/Pe)^(2/ke)-(Pc/Pe)^((ke+1)/ke)));
                    dme=Ae*Ce*rhoe*Fe/omega;
                    C2=Te;
                end % end exhaust Pe/Pc check
                if Pi/Pc>=Cri %inward choked flow
                    Fi=(2/(ki+1))^((ki+1)/2/(ki-1));
                    dmi=Ai*Ci*rhoi*Fi/omega; %intake valve mass flow rate dm/dtheta (kg/deg)
                    C1=Ti;
                elseif Pi/Pc<Cri&&Pi/Pc>=1 %inward nonchoked flow
                    Fi=sqrt(2/(ki-1)*((Pc/Pi)^(2/ki)-(Pc/Pi)^((ki+1)/ki)));
                    dmi=Ai*Ci*rhoi*Fi/omega;
                    C1=Ti;
                elseif Pc/Pi>=Cri %outward choked flow
                    Fi=-(2/(ki+1))^((ki+1)/2/(ki-1));
                    dmi=Ai*Cc*rhoc*Fi/omega; %
                    C1=Tc;
                else %outward nonchoked flow
                    Fi=-sqrt(2/(ki-1)*((Pi/Pc)^(2/ki)-(Pi/Pc)^((ki+1)/ki)));
                    dmi=Ai*Cc*rhoc*Fi/omega;
                    C1=Tc;
                end % end intake Pi/Pc check
            end % end overlap flow calc
            % only intake valve open
            if theta>EVC && theta<=IVC
                if Pi/Pc>=Cri %inward choked flow
                    Fi=(2/(ki+1))^((ki+1)/2/(ki-1));
                    dmi=Ai*Ci*rhoi*Fi/omega; %intake valve mass flow rate dm/dtheta (Kg/deg)
                    C1=Ti;
                    C2=0;
                elseif Pi/Pc<Cri&&Pi/Pc>=1 %inward nonchoked flow
                    Fi=sqrt(2/(ki-1)*((Pc/Pi)^(2/ki)-(Pc/Pi)^((ki+1)/ki)));
                    dmi=Ai*Ci*rhoi*Fi/omega;
                    C1=Ti;
                    C2=0;
                elseif Pc/Pi>=Cri %outward choked flow
                    Fi=-(2/(ki+1))^((ki+1)/2/(ki-1));
                    dmi=Ai*Cc*rhoc*Fi/omega; %
                    C1=Tc;
                    C2=0;
                else %outward nonchoked flow
                    Fi=-sqrt(2/(ki-1)*((Pi/Pc)^(2/ki)-(Pi/Pc)^((ki+1)/ki)));
                    dmi=Ai*Cc*rhoc*Fi/omega;
                    C1=Tc;
                    C2=0;
                end % end Pi/Pc check
                dme=0; % no mass flow through exhaust valve
            end % end intake flow calc
            term1=-ka*Pc*dvol/vol;
            term2= Rs*ka*C1*dmi/vol;
            term3= Rs*ka*C2*dme/vol;
            yprime(1,1)=term1+term2+term3; %dP/dtheta
            yprime(2,1)=dmi+dme; %net mass inflow
            yprime(3,1)=Pc*dvol; %PdV/dtheta work
        end %end of function rates
    end % end of function integrate
end % end of main program function intake exhaust