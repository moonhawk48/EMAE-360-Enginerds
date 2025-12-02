

clc
clear all
close all

corsa =0.108;          %Stroke m
ales = 0.0864;          %Bore m
lbiel =   0.21;           %Connecting rod length m
rapcom =   10;             %Compression ratio %
dporta =  0.05;           %Intake valve diameter m
anta =     30;             %Intake valve opening advance ° 
rita =  20;             %Intake valve closing delay °
dports = 0.04;          %Exhaust valve diameter m
ants =  60;             %Exhaust valve opening advance ° 
rits =   25;             %Exhaust valve closing delay °
nc =     2;              %Number of engine cylinders []

% Si introducono le condizioni operative del motore 


alfa =  14.5 ;            %Air-fuel ratio ma/mc
rpm =   800 ;             %Rotation speed rpm
thant =  -20 ;            %Ignition advance °
dthb =  50 ;              %Angular combustion duration °
pasp = 0.95 ;            %Intake manifold pressure bar
tasp = 300 ;             %Intake manifold temperature K       
psca =   1.05 ;            %Exhaust manifold pressure bar
tsca =    600 ;             %Exhaust manifold temperature Exhaust temperature K
twall =  450 ;            %Temperature of the heat exchange surfaces between the mixture and the thermal unit K
 
 % Tipo di motore
 
 id = 0;            %Set id = 1 if you want to calculate the cycle for a diesel engine, otherwise it will be acc_com
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Fine parametri di controllo                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Si riportano delle relazioni notevoli valide per qualsiasi motore e che
% risulteranno utili nei calcoli successivi

 apisto = pi * ales^2 / 4 ;              %Piston head surface area m^2
 rcl = ales / 2 ;                        %Cylinder radius m
 rman = corsa / 2 ;                      %Crank radius m
 vvs = apisto * corsa ;                  %Unit displacement m^3
 rsul = rman / 2 / lbiel ;               %mu means (parameter related to the crank mechanism) /
 vvc = vvs / (rapcom - 1) ;              %Combustion chamber volume m^3
 fgr = pi/180;                           %Conversion of degrees to radians
 



up = rpm * corsa / 30;

%Evaluate the valve lift laws

 durataa = 180 + anta + rita ;          %Angular duration of the intake °
 sezionea = pi * dporta^2 / 4 ;         %Reference section of the intake valves m^2
 
[Tasp,Aasp]= OdLibv2.Vasp(durataa,sezionea,anta);

duratas = 180 + ants + rits ;           %Angular duration of the exhaust °
seziones = pi * dports^2 / 4 ;          %Reference section of the exhaust valves m^2

[Tsca,Asca]= OdLibv2.Vsca(duratas,seziones,ants);

Taspsca = [Tasp' Tsca'];                %Lift angle matrix
Aaspsca = [Aasp' Asca'];                %Angle section matrix Lifted

% Disegno le leggi calcolate

figure (1)

plot  (Tsca,Asca)
hold on
plot  (Tasp,Aasp)
xlabel('Theta')
ylabel('Valve lift')

% Caratteristiche della benzina

PCIb = 10500 * 4.186 ;                   %Gasoline calorific value kJ/kg
k = 1.35 ;                               %Polytropic exponent /
r = 297 ;                                %Universal ideal gas constant kJ/kgK
cp = k * r / (k - 1) / 1000 ;            %Specific heat at constant pressure kJ/kgK
cv = r / (k - 1) / 1000  ;               %Specific heat at constant volume kJ/kgK
gam = cp/cv;                             %Gamma ratio
aa = 5 ;                                 %Wiebe model constants
e = 2 ;                                  

% Condizioni ambiente esterno

pamb = 101325;                           %Ambient pressure Pa
tamb = 300;                              %Ambient temperature K
hamb = cp * tamb;                        %Ambient specific enthalpy kJ/kg
roamb = pamb / (r * tamb);               %Ambient density kg/m^3
 
% Condizioni aspirazione                

pasp = pasp * 100000;                   %Conversion pasp to Pa
hasp = cp * tasp;                       %Intake enthalpy
roasp = pasp / (r * tasp);              %Intake density

% Condizione scarico
      
psca = psca * 100000;                   %Conversion psca to Pa
hsca = cp * tsca;                       %Exhaust enthalpy
rosca = psca / (r * tsca);              %Exhaust density

% Struttura ciclica del modello

 %ncicli = 3;                             %Number of cycles to be performed
 tstart = -180 + rita;                   %Start angle of calculation
 tstop = tstart + 720;                   %End angle Cycle calculation
 %tstopt = tstart + 720 * ncicli;         %End calculation angle
 %idciclo = 1;                            %Cycle index
 dteta = 0.2;                            %Angular calculation step
 neq = 3;                                %Number of equations to integrate
 %tetacicli = tstart:dteta:tstopt;        %Vector containing all the angles analyzed in the cycles
 tetaciclo = tstart:dteta:(tstart+720);  %Vector containing all the angles analyzed in the cycle
 siz = size (tetaciclo,2);               %Number of crank angles analyzed

% Condizioni iniziali

etab = 0.95;                             %Combustion efficiency
etam = 0.9;                              %Mechanical efficiency
pnot = pasp;                             %Known pressure
tnot = tasp;                             %Known temperature
hnot = cp * tnot;                        %Known enthalpy
unot = cv * tnot;                        %Known internal energy
teta = tstart;                           %Start calculation angle
vnot = OdLibv2.Vdv (fgr,teta,vvc,vvs,rsul);        %Calls the library function to calculate the available fluid volume.
mnot = pnot * vnot / r / tnot;           %Known Mass
mref = pamb * vvs / r / tamb;            %Reference Mass (referred to displacement)

mb = 0 ;                                 %Burn Mass

Odxb = zeros(siz,1);          %Initialization of the xb vector
x = [mb tnot mnot];           %Initialization of the parameter vector to be calculated using Merson integration
Ox = ones (siz,3);            %Initialization of the parameter vector matrix to be calculated using Merson integration
Ox (1,:)=[x];                 


for i = 1:siz-1                                  %First Cycle
    teta = tetaciclo(i);                         %Indexing of Crank Angle
    x = Ox(i,:);                                 %Selection of Starting Values
    
    %Call the Merson function from the library
    [x] = OdLibv2.Merson (neq,teta,dteta,x,fgr,vvc,vvs,rsul,apisto,k,r,pnot,vnot,twall,rpm,thant,dthb,id,alfa,PCIb,cv,cp,pasp,psca,hasp,hsca,roasp,rosca,Taspsca,Aaspsca,ales,up,tnot,aa,e);
    Ox(i+1,:)= [x];                             %Result of Merson interaction placed in row i+1
    
    %Call the Rates function from the library to calculate the rate of
    %heat release at indexed theta
    D = teta+dteta;                             
    dxb = OdLibv2.rates (D,x,fgr,vvc,vvs,rsul,apisto,k,r,pnot,vnot,twall,rpm,thant,dthb,id,alfa,PCIb,cv,cp,pasp,psca,hasp,hsca,roasp,rosca,Taspsca,Aaspsca,ales,up,tnot,aa,e);
   Odxb(i+1,1) = dxb(1);    
end

Ox (1,:)=[x];                                   %Riinizializza la matrice dei parametri calcolati assegnando come primo vettore riga l'ultimo vettore del ciclo precedente
Ox (1,1) = 0;                                   %Riporta a zero la massa bruciata
Odxb = zeros(siz,1);
for i = 1:siz-1                                  %Second Cycle
    teta = tetaciclo(i);
    x = Ox(i,:);
    [x] = OdLibv2.Merson (neq,teta,dteta,x,fgr,vvc,vvs,rsul,apisto,k,r,pnot,vnot,twall,rpm,thant,dthb,id,alfa,PCIb,cv,cp,pasp,psca,hasp,hsca,roasp,rosca,Taspsca,Aaspsca,ales,up,tnot,aa,e);
    Ox(i+1,:)= [x];
    D = teta+dteta;
    dxb = OdLibv2.rates (D,x,fgr,vvc,vvs,rsul,apisto,k,r,pnot,vnot,twall,rpm,thant,dthb,id,alfa,PCIb,cv,cp,pasp,psca,hasp,hsca,roasp,rosca,Taspsca,Aaspsca,ales,up,tnot,aa,e);
   Odxb(i+1,1) = dxb(1);
end

Ox (1,:)=[x];
Ox (1,1) = 0;
Odxb = zeros(siz,1);
for i = 1:siz-1                                  %Third Cycle
    teta = tetaciclo(i);
    x = Ox(i,:);
    [x] = OdLibv2.Merson (neq,teta,dteta,x,fgr,vvc,vvs,rsul,apisto,k,r,pnot,vnot,twall,rpm,thant,dthb,id,alfa,PCIb,cv,cp,pasp,psca,hasp,hsca,roasp,rosca,Taspsca,Aaspsca,ales,up,tnot,aa,e);
    Ox(i+1,:)= [x];
    D = teta+dteta;
    dxb = OdLibv2.rates (D,x,fgr,vvc,vvs,rsul,apisto,k,r,pnot,vnot,twall,rpm,thant,dthb,id,alfa,PCIb,cv,cp,pasp,psca,hasp,hsca,roasp,rosca,Taspsca,Aaspsca,ales,up,tnot,aa,e);
   Odxb(i+1,1) = dxb(1);
end


for i = 1:siz                                 
    teta = tetaciclo(i);
    %Call the crank law from the library to find V and dv
    %as theta varies
[V(i),dv(i)] = OdLibv2.Vdv (fgr,teta,vvc,vvs,rsul);
%Use mass, temperature, and volume to calculate pressure
%using the ideal gas law
 P (i)= (r*Ox(i,3)*Ox(i,2)/V(i))/100000;

end

    % Ciclo di pompaggio
    
    Pdue = [P P];
    Vdue = [V V];
    
    
% Valuto il lavoro utile indicato 

Lui = zeros (siz,1);                            %Initialization

for i = 1:siz-1
    
Lui(i+1) = Lui(i) + P(i)*100000 * dv(i) * dteta;

end


Lui = Lui/1000;                                     %Indicated useful work

% Calcolo i dati tecnici del motore
    
Lue = Lui * etam;                                   %Effective work
display 'Effective Work'
potkw = nc * Lue(siz) * rpm / 60 / 2                %Power in kW
display 'Torque'
copnm = potkw * 1000 * 60 / 2 / pi / rpm            %Torque     
display 'Avg Pressure'
pmi = Lui(siz) / vvs /100                           %Average pressure
mtrap = Ox (siz,3);
mf = mtrap / (alfa + 1);
display 'Specific fuel consumption'
cs = (mf * 1000) / (Lue(siz) / 3600) / etab         %Specific consumption 
display 'Riempimento'
lambdav = mtrap / mref                              %Filling 



 Odxb= Odxb*100000/4;                                            %Scale adjustment
  
  % Salva i dati su file esterno

%   save Calcolo.mat;
   
  % plot
  figure (5)                                                     %Heat release law
  plot (tetaciclo(250:1000),Ox(250:1000,1)/Ox(250:1000,3)); 
  
 hold on                                                        %Release rate
  plot (tetaciclo(250:1000),Odxb(250:1000,1),'DisplayName','Release Rate');
   xlabel('Theta')
   legend()

  
  
 figure (7)                                             %P-V diagram
  plot(V,P); 
        xlabel('Volume')
        ylabel('Pressure')

   figure (8)                                                   %Progression of indicated work during the cycle
  plot (tetaciclo,Lui);
      xlabel('Theta')
  ylabel('Indicated Work')
  
   figure (2)                                 %P-theta cycle
  plot (tetaciclo,P);  
    xlabel('Theta')
  ylabel('Pressure')
%     
%     subplot (3,1,3)
  figure (3)                                %Volume variation law

    plot (tetaciclo,V)
          xlabel('Theta')
        ylabel('Volume')
    
    % Pumping Cycle
    
    Pdue = [P P];
    Vdue = [V V];
    
    figure (10)
    plot ( Vdue (1700:3800),Pdue(1700:3800)  )
title('Pumping Cycle')
    
     figure (4)                                         %T-theta graph
  plot (tetaciclo,Ox(:,2));
  xlabel('Theta')
  ylabel('Temperature')
  
 figure (6)                                         %Cylinder involute mass graph
  plot (tetaciclo,Ox(:,3));
    xlabel('Theta')
  ylabel('Cylinder Unburned Mass')

  % Added by Sean for flywheel calcs (need some torque stuff)

  % Instant power at crank angle and instantaneous torque at crank angle
  for z=1:siz
    powKW(z) = nc * Lue(z) * rpm / 60 / 2;
    tqNM(z) = powKW(z) * 1000 * 60 / 2 / pi / rpm;
end

% Mean torque
meanTQ = mean(tqNM)

% Subtract mean torque from torque
modTQ = tqNM - meanTQ;

% Delta E for flywheel is integral of T - T mean
DeltaE = trapz(modTQ)

    
