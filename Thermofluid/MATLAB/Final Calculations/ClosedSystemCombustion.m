function [out] = ClosedSystemCombustion()
%Calculates the closed system portion of the otto cycle
clc
workspace;

% Engine geometry and operating conditions
R = 10; % Compression ratio - R
B = .108; % Bore - B (m)
S = .0864; % Stroke - S (m)
EPS = 0.25; % Half stroke to rod ratio - EPS
RPM = 6000; % Engine speed - RPM
BLOWBY = 0.8; % Blowby coefficient - Dependent on piston ring setup
HEAT = 200; %Initial ht coefficient
PHI = 1; % Equivalence ratio - PHI
F = 0.05; % Residual fraction - F
TW = 500; % Wall temperature - TW
fuel_type = 2; % gasoline
FS = 0.06548; % stoichiometric fuel-air ratio for gasoline
A0 = 47870; %
%inlet temperature and pressure
T1 = 350;
P1 = 101.172; % kPa


to_ppm = 10^6; % convert from mass fraction to ppm
MW_NO = 30;  % molecular weight of NO, g/mol

THETAS = -25; % Start of heat release (deg ATDC)
THETAB = 30; % Burn angle (deg)
IVC = -180;
EVO = 180;


%inputs: Input temp, input pressure, RPM, Wall temp, charge mass, IVC, EVO

    % R= varargin{1};
    % B= varargin{2};
    % S= varargin{3};
    %T1 = varargin{1};
    % P1 = varargin{5};
    %RPM = varargin{2};
    % TW = varargin{7};
    % MNOT = varargin{8};
    % IVC = varargin{9}-360;
    % EVO = varargin{10}-360;

%engine geometries    
OMEGA = RPM*pi/30; %RPM to rad
THETA = IVC; %IVC
THETA_END = EVO; % EVO
DTHETA = 1; %step size
THETAE = THETA+DTHETA;
finalStep = floor((THETA_END+abs(THETA))/10);
Up = S*OMEGA/pi; %mean piston speed

%engine geometry function call
[ VOL, X, EM ] = auxiliary( THETA );
NNOX = THETAB/DTHETA;
NY = 6+NNOX;
Y = zeros(NY,1);
Y(1) = P1;
Y(2) = nan;
Y(3) = T1;
%initial cylinder energy
[~, ~, ~, ~, vU, ~, ~, ~, ~, ~] = farg( Y(3), Y(1), PHI, F, fuel_type );
MNOT = VOL/vU;
M = EM*MNOT;

%zero arrays for outputs
NN = 36*10;
SAVE.THETA = zeros( NN, 1 );
SAVE.VOL = zeros( NN, 1 );
SAVE.T = zeros(NN, 1 );
SAVE.P = zeros( NN, 1 );
SAVE.MDOTFI = zeros( NN, 1 );
SAVE.NOX = zeros(NN,5);
SAVE.NOe = zeros(NN,1);
%display outputs
fprintf( 'THETA    VOLUME   BURN FRAC  PRESSURE    BURN TEMP   UNBURNED T     WORK    HEAT LOSS    MASS      H-LEAK      NOx\n' );
fprintf( ' deg      cm^3       --       kPa           K              K         J           J          g          J       ppm\n' );

fprintf('%7.1f   %6.1f   %3.3f      %6.1f       %6.1f         %6.1f     %5.0f     %5.0f       %5.3f    %5.2f     %6.1f\n', ...
    THETA, VOL*1000000, X, Y(1), Y(2), Y(3), Y(4)*1000, Y(5)*1000, M*1000, Y(6)*1000, 0.0 );
II = 1;
%loop through combustion and expansion angles
for III=1:finalStep      
    for JJJ=1:10
        NOe_save = 0;
        %perform integration of differentials
        [ Y ] = integrate( THETA, THETAE, Y );
        [ VOL, X, EM ] = auxiliary( THETA );
        %calculate ht coefficient
        HEAT = Woschni();

        %mass inducted
        M = EM*MNOT;
        THETA=THETAE;
        THETAE=THETA+DTHETA;
        
        % save data for plotting later
        SAVE.THETA(II) = THETA;
        SAVE.VOL(II) = VOL;
        SAVE.P(II) = Y(1);
        SAVE.TB(II) = Y(2);
        SAVE.TU(II) = Y(3);
        SAVE.W(II) = Y(4);
        SAVE.Q(II) = Y(5);
        SAVE.X(II) = X;

        SAVE.NOX(II,:) = [ Y(6+1), Y(round(6+0.25*NNOX)),...
            Y(round(6+0.5*NNOX)), Y(round(6+0.75*NNOX)), Y(6+NNOX) ]*to_ppm;
        SAVE.NOe(II) = NOe_save*to_ppm; % save equilibrium NOx mass fraction for plotting
        
        II=II+1;

        if ( THETAS >= THETA && THETAS < THETAE )
           Y(2) = tinitial( Y(1), Y(3), PHI, F );
        end
        
        if ( THETA > THETAS + THETAB )
            Y(3) = nan;

        end 
        fprintf('%7.1f   %6.1f   %3.3f      %6.1f       %6.1f         %6.1f     %5.0f     %5.0f       %5.3f    %5.2f     %6.1f\n', ...
        THETA, VOL*1000000, X, Y(1), Y(2), Y(3), Y(4)*1000, Y(5)*1000, M*1000, Y(6)*1000, Y(7)*to_ppm );
    
    end
     
end

% integrate total NOx value
NOX_ppm = 0;
for nn=1:NNOX;
    THETA = THETAS + (nn-1)/(NNOX-1)*THETAB;
    dxbdtheta = 0.5*sin(pi*(THETA-THETAS)/THETAB)*pi/THETAB;
    dxb = dxbdtheta*DTHETA;
    NOX_ppm = NOX_ppm + Y(6+nn)*dxb*to_ppm;
end
%calculate thermal efficiency and IMEP
ETA = Y(4)/MNOT*(1+PHI*FS*(1-F))/PHI/FS/(1-F)/A0;
IMEP = Y(4)/(pi/4*B^2*S);
fprintf('ETA = %1.4f  IMEP = %7.3f kPa  NOx = %6.1f ppm\n', ETA, IMEP, NOX_ppm );
% T_out = mean([SAVE.TU,SAVE.TB],2,'omitnan')
% P_out = SAVE.P;
W_total = Y(4);
Q_total = Y(5);
out = [W_total,Q_total];
%if calling outside of excel, create plots
if ( nargin == 2 )

    sTitle = sprintf('Homogenous 2 zone, methane, PHI=%.2f F=%.2f RPM=%.1f\nETA=%.3f  IMEP=%.2f kPa  NOx=%.1f ppm ', PHI, F, RPM, ETA, IMEP, NOX_ppm );

    figure;
    plot( SAVE.THETA, SAVE.X, 'linewidth',2 );
    set(gca,'fontsize',18,'linewidth',2,'Xlim',[-100 100]); 
    xlabel( '\theta','fontsize',18);
    ylabel('burn fraction','fontsize',18);
    print -deps2 NOxfrac
    title( sTitle );
     
    figure;
    plot(  SAVE.THETA, SAVE.P,'linewidth',2 );
    set(gca,'fontsize',18,'linewidth',2,'Xlim',[-100 100]); 
    xlabel( '\theta','fontsize',18);
    ylabel('pressure (kPa)','fontsize',18);
    print -deps2 NOxpress
    
    figure;
   plot( SAVE.THETA, SAVE.TU, '-',SAVE.THETA, SAVE.TB,'--','linewidth',2 );
    set(gca,'fontsize',18,'linewidth',2,'Xlim',[-100 100]); 
    xlabel( '\theta','fontsize',18);
    ylabel( 'temperature (K)', 'fontsize',18); 
    legend('Unburned','Burned', 'Location', 'SouthEast');
   

    figure;
    plot( SAVE.THETA, SAVE.NOX, SAVE.THETA, SAVE.NOe,'linewidth',2  );
     set(gca,'fontsize',18,'linewidth',2,'Xlim',[-100 100]); 
    xlabel('\theta','fontsize',18);
    ylabel('NOx (ppm)','fontsize',18);
    axis( [ THETAS, 110, 0, max(SAVE.NOe)*1.1 ] );
     print -deps2 NOxNOx
    legend( 'X=0', 'X=0.25', 'X=0.5', 'X=0.75', 'X=1', 'Equilibrium', 'Location', 'NorthEast' );
    %title( sTitle );
    
    disp(Up)

end

%helper function for burn temperature
function [ TB ] = tinitial( P, TU, PHI, F )
    TB = 2000;
    [~, HU,~, ~, ~, ~, ~, ~, ~, ~] = farg( TU, P, PHI, F, fuel_type );
    for ITER=1:50
        [ierr, ~, HB,~, ~, ~, ~, CP, ~, ~, ~] = ecp( TB, P, PHI, fuel_type );
        if ( ierr ~= 0 )
            fprintf('Error in ECP(%g, %g, %g): %d\n', TB, P, PHI, ierr )
        end
        DELT = +(HU-HB)/CP;
        TB = TB + DELT;
        if ( abs(DELT/TB) < 0.001 )
            break;
        end
    end
end

%helper function for engine kinematics
function [ VOL, X, EM ] = auxiliary( THETA )
    VTDC = pi/4*B^2*S/(R-1); % m3
    VOL = VTDC*(1 + (R-1)/2*(1-cosd(THETA) + 1/EPS*(1-sqrt(1-(EPS*sind(THETA))^2))));
    X = 0.5*(1-cos(pi*(THETA-THETAS)/THETAB));
    if ( THETA <= THETAS )
        X = 0.;
    end
    if ( THETA >= THETAS+THETAB )
        X = 1.; 
    end
    EM = exp(-BLOWBY*(THETA*pi/180 + pi)/OMEGA);
end

%integration of differential functions
function [Y] = integrate( THETA, THETAE, Y )

[TT, YY ] = ode23( @rates, [ THETA, THETAE ], Y );
for J=1:NY,
    Y(J) = YY(length(TT),J);
end

    function [ YPRIME ] = rates( THETA, Y )
        YPRIME = zeros(NY,1);
        
        [ VOL, X, EM ] = auxiliary( THETA );
        M = EM*MNOT;
        DUMB = sqrt(1-(EPS*sind(THETA))^2);
        DV = pi/8*B^2*S*sind(THETA)*(1+EPS*cosd(THETA)/DUMB);
        AA = (DV + VOL*BLOWBY/OMEGA)/M;
        C1 = HEAT*(pi*B^2/2 + 4*VOL/B)/OMEGA/M/1000;
        C0 = sqrt(X);
        
        P = Y(1);
        TB = Y(2);
        TU = Y(3);
        
        % three different computations are required depending upon the size
        % of the mass fraction burned
        if ( X > 0.999 )
            %  EXPANSION
            [ierr, YB, HL, ~, ~, VB, ~, CP, ~, DVDT, DVDP] = ecp( TB, P, PHI, fuel_type );
            if ( ierr ~= 0 )
                fprintf('Error in ECP(%g, %g, %g): %d\n', TB, P, PHI, ierr );
            end                       
            
            BB = C1/CP*DVDT*TB*(1-TW/TB);
            CC = 0;
            DD = 1/CP*TB*DVDT^2 + DVDP;
            EE = 0;
            
            YPRIME(1) = (AA + BB + CC)/(DD + EE);
            YPRIME(2) = -C1/CP*(TB-TW) + 1/CP*DVDT*TB*YPRIME(1);
            YPRIME(3) = 0;
            
        elseif ( X > 0.001 )
            %  COMBUSTION
            [~, HU, ~, ~, VU, ~, CPU, ~, DVDTU, DVDPU] = farg( TU, P, PHI, F, fuel_type );
            [ierr, YB, HB, ~, ~, VB, ~, CPB, ~, DVDTB, DVDPB] = ecp( TB, P, PHI, fuel_type );
            if ( ierr ~= 0 )
                fprintf('Error in ECP(%g, %g, %g): %d\n', TB, P, PHI, ierr );
            end
            
            BB = C1*(1/CPB*TB*DVDTB*C0*(1-TW/TB) + 1/CPU*TU*DVDTU*(1-C0)*(1-TW/TU));
            DX = 0.5*sin( pi*(THETA-THETAS)/THETAB )*180/THETAB;
            CC = -(VB-VU)*DX - DVDTB*(HU-HB)/CPB*(DX-(X-X^2)*BLOWBY/OMEGA);
            DD = X*(VB^2/CPB/TB*(TB/VB*DVDTB)^2 + DVDPB);
            EE = (1-X)*(1/CPU*TU*DVDTU^2 + DVDPU);
            HL = (1-X^2)*HU + X^2*HB;
            
            YPRIME(1) = (AA + BB + CC)/(DD + EE);
            YPRIME(2) = -C1/CPB/C0*(TB-TW) + 1/CPB*TB*DVDTB*YPRIME(1) + (HU-HB)/CPB*(DX/X - (1-X)*BLOWBY/OMEGA);
            YPRIME(3) = -C1/CPU/(1+C0)*(TU-TW) + 1/CPU*TU*DVDTU*YPRIME(1);
            
        else
            %  COMPRESSION      
            [~, HL, ~, ~, ~, ~, CP, ~, DVDT, DVDP] = farg( TU, P, PHI, F, fuel_type );
            
            BB = C1*1/CP*TU*DVDT*(1-TW/Y(3));
            CC = 0;
            DD = 0;
            EE = 1/CP*TU*DVDT^2 + DVDP;
            
            YPRIME(1) = ( AA + BB + CC )/(DD + EE);
            YPRIME(2) = 0;
            YPRIME(3) = -C1/CP*(Y(3)-TW) + 1/CP*Y(3)*DVDT*YPRIME(1);
        end
        
        % save equilibrium NO concentration
        NOe_save = 0;
        if ( X > 0.001 )
            % save the equilibrium NO concentration as mass fraction
            N_V = (1000*P)/(8.314*TB)*(1/100)^3; % calculate total molar concentration [mol/cm^3]
            NOe_save = YB(10)*N_V*MW_NO*VB*1000; % NO mass fraction
        end
        % common to all cases
        YPRIME(4) = Y(1)*DV;
        
        YPRIME(5) = 0;
        if ( ~isnan(TB) )
            YPRIME(5) = YPRIME(5) + C1*M*C0*(TB-TW);
        end
        if ( ~isnan(TU) )
            YPRIME(5) = YPRIME(5) + C1*M*(1-C0)*(TU-TW);
        end
        
        YPRIME(6) = BLOWBY*M/OMEGA*HL;

        % 1/omega is s/rad, so convert to s/deg
        for JJ=1:NY,
            YPRIME(JJ) = YPRIME(JJ)*pi/180;
        end
    end
end

function [h] = Woschni()
        P_ref = (P1)*(VOL(1)/VOL)^1.3;
        if ( THETA > THETAS + THETAB )
           T_ht = Y(2);
        else
          T_ht = Y(3);
        end
        
        if (THETA < THETAS) 
            U = 2.28*Up; 
        else
            term4 = T1*(VOL/VOL(1))*(Y(1)-P_ref)/P1;
            U = 2.28*Up + 0.00324*term4; % Woschni vel (m/s);
        end
        h = 3.26*((Y(1))^0.8)*(U^0.8)*(B^-0.2)*T_ht^-0.55;
end
end
