                                                                          
  %%%%%%%%%%%%%%    %%%%%%%%%%%      %%%     %%%      %%%%%%%%%%%         
       %%%          %%%              %%%%   %%%%      %%%                 
       %%%          %%%              %%% % % %%%      %%%                 
       %%%          %%%              %%%  %  %%%      %%%%%%%%%%%         
       %%%          %%%              %%%     %%%              %%%         
       %%%          %%%              %%%     %%%              %%%         
       %%%rue       %%%%%%%%%%%ycle  %%%     %%%otor  %%%%%%%%%%%imulator 
                                                                    
       % Mirko Casale, Marco Marzocchi, Dario Santoro, Simone Talento %
       
       % Conversione in Matlab di cicli0D (VB) del prof. Fabio Bozza %
       
                                    %2019%
                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Questo script permette di calcolare e rappresentare il ciclo reale di un%
% motore termico al variare dei parametri evidenziati di seguito.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script allows you to calculate and represent the real cycle of a   %
% ICE as the parameters shown below vary.                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    

% Inizializzazione delle spazio di calcolo

clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Inizio dei parametri di controllo                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Si introducono le caratteristiche geometriche del motore


corsa =0.072;          %Corsa m
ales = 0.087;          %Alesaggio m
lbiel =   0.15;           %Lunghezza biella m
rapcom =   11;             %Rapporto di compressione %
dporta =  0.03;           %Diametro valvola di aspirazione m
anta =     30;             %Anticipo apertura valvola di aspirazione ° 
rita =  30;             %Ritardo alla chiusura della valvola di aspirazione °
dports = 0.025;          %Diametro valvola di scarico m
ants =  60;             %Anticipo apertura valvola di scarico ° 
rits =   30;             %Ritardo alla chiusura della valvola di scarico °
nc =     4;              %Numero di cilindri del motore []

% Si introducono le condizioni operative del motore 

 

alfa =  14.5 ;            %Rapporto aria combustibile ma/mc
rpm =   3000 ;             %Velocità di rotazione giri/min
thant =  -20 ;            %Anticipo all'accensione °
dthb =  50 ;              %Durata angolare della combustione °
pasp = 0.95 ;            %Pressione al collettore di aspirazione bar
tasp = 300 ;             %Temperatura al collettore di aspirazione K       
psca =   1.05 ;            %Pressione al collettore di scarico bar
tsca =    600 ;             %Temperatura al collettore di scarico K
twall =  450 ;            %Temperatura delle superfici di scabio termico tra miscela e gruppo termico K   
 
 % Tipo di motore
 


 
 id = 0;            %Impostare id = 1 se si vuole calcolare il ciclo per un motore disel altrimenti sarà acc_com
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Fine parametri di controllo                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Si riportano delle relazioni notevoli valide per qualsiasi motore e che
% risulteranno utili nei calcoli successivi

 apisto = pi * ales^2 / 4 ;              %Superficie testa del pistone m^2
 rcl = ales / 2 ;                        %Raggio del cilindro m
 rman = corsa / 2 ;                      %Raggio della manovella m
 vvs = apisto * corsa ;                  %Cilindrata unitaria m^3
 rsul = rman / 2 / lbiel ;               %mu mezzi (paramtro legato al manovellismo) /
 vvc = vvs / (rapcom - 1) ;              %Volume camera di combustione m^3
 fgr = pi/180;                           %Conversione gradi in radianti 
 



up = rpm * corsa / 30;

%Valuto le leggi di alzata delle valvole

 durataa = 180 + anta + rita ;          %Durata angolare dell'aspirazione °
 sezionea = pi * dporta^2 / 4 ;         %Sezione di riferimento valvole di aspirazione m^2
 
[Tasp,Aasp]= OdLibv2.Vasp(durataa,sezionea,anta);

duratas = 180 + ants + rits ;           %Durata angolare dello scarico °
seziones = pi * dports^2 / 4 ;          %Sezione di riferimento valvole di scarico m^2

[Tsca,Asca]= OdLibv2.Vsca(duratas,seziones,ants);

Taspsca = [Tasp' Tsca'];                %Matrice degli angoli di alzata
Aaspsca = [Aasp' Asca'];                %Matrice delle sezioni di alzata

% Disegno le leggi calcolate

figure (1)

plot  (Tsca,Asca)
hold on
plot  (Tasp,Aasp)

% Caratteristiche della benzina

PCIb = 10500 * 4.186 ;                   %Potere calorifico benzina kJ/kg
k = 1.35 ;                               %Esponente della politropica /
r = 297 ;                                %Costante universale dei gas perfetti kJ/kgK
cp = k * r / (k - 1) / 1000 ;            %Calore specifico a pressione costante kJ/kgK 
cv = r / (k - 1) / 1000  ;               %Calore specifico a volume costante kJ/kgK 
gam = cp/cv;                             %Rapporto gamma
aa = 5 ;                                 %Costanti del modello Wiebe
e = 2 ;                                  

% Condizioni ambiente esterno

pamb = 101325;                           %Pressione ambiente Pa      
tamb = 300;                              %temperatura ambiente K           
hamb = cp * tamb;                        %Entalpia specifica ambiente kJ/kg
roamb = pamb / (r * tamb);               %Densità ambiente kg/m^3
 
% Condizioni aspirazione                

pasp = pasp * 100000;                   %Conversione pasp in Pa
hasp = cp * tasp;                       %Entalpia aspirazione
roasp = pasp / (r * tasp);              %Densità aspirazione

% Condizione scarico
      
psca = psca * 100000;                   %Conversione psca in Pa
hsca = cp * tsca;                       %Entalpia scarico
rosca = psca / (r * tsca);              %Densità scarico

% Struttura ciclica del modello

 %ncicli = 3;                             %Numero di cicli da eseguire
 tstart = -180 + rita;                   %Angolo inizio calcolo
 tstop = tstart + 720;                   %Angolo fine calcolo ciclo
 %tstopt = tstart + 720 * ncicli;         %Angolo fine calcolo 
 %idciclo = 1;                            %Indice del ciclo
 dteta = 0.2;                            %Step angolare di calcolo
 neq = 3;                                %Numero di equazioni da integrare
 %tetacicli = tstart:dteta:tstopt;        %Vettore che contiene tutti gli angoli analizzati nei cicli
 tetaciclo = tstart:dteta:(tstart+720);  %Vettore che contiene tutti gli angoli analizzati nel ciclo
 siz = size (tetaciclo,2);               %Numero degli angoli di manovella analizzati                  

% Condizioni iniziali

etab = 0.95;                             %Rendimento di combustione
etam = 0.9;                              %Rendimento meccanico
pnot = pasp;                             %Pressione nota
tnot = tasp;                             %Temperatura nota
hnot = cp * tnot;                        %Entalpia nota
unot = cv * tnot;                        %Energia interna nota
teta = tstart;                           %Angolo di inizio calcolo
vnot = OdLibv2.Vdv (fgr,teta,vvc,vvs,rsul);        %Richiama la funzione della libreria per valutare il volume a disposizione del fluido.
mnot = pnot * vnot / r / tnot;           %Massa nota
mref = pamb * vvs / r / tamb;            %Massa di riferimento (riferita alla cilindrata)

mb = 0 ;                                 %Massa bruciata

Odxb = zeros(siz,1);          %Inizializzazione del vettore xb
x = [mb tnot mnot];           %Inizializzazione del vettore dei parametri da calcolare tramite integrazione Merson
Ox = ones (siz,3);            %Inizializzazione della matrice dei vettori dei parametri da calcolare tramite integrazione Merson
Ox (1,:)=[x];                 


for i = 1:siz-1                                  %Primo ciclo
    teta = tetaciclo(i);                         %Indicizzazione angolo di manovella
    x = Ox(i,:);                                 %Selezione dei valori di partenza
    
    %Richiama la funzione Merson dalla libreria
    [x] = OdLibv2.Merson (neq,teta,dteta,x,fgr,vvc,vvs,rsul,apisto,k,r,pnot,vnot,twall,rpm,thant,dthb,id,alfa,PCIb,cv,cp,pasp,psca,hasp,hsca,roasp,rosca,Taspsca,Aaspsca,ales,up,tnot,aa,e);
    Ox(i+1,:)= [x];                             %Risultato interazione Merson collocato nella riga i+1
    
    %Richiama la funzione Rates dalla libreria per calcolare la velocità di
    %rilascio del calore al teta indicizzato
    D = teta+dteta;                             
    dxb = OdLibv2.rates (D,x,fgr,vvc,vvs,rsul,apisto,k,r,pnot,vnot,twall,rpm,thant,dthb,id,alfa,PCIb,cv,cp,pasp,psca,hasp,hsca,roasp,rosca,Taspsca,Aaspsca,ales,up,tnot,aa,e);
   Odxb(i+1,1) = dxb(1);    
end

Ox (1,:)=[x];                                   %Riinizializza la matrice dei parametri calcolati assegnando come primo vettore riga l'ultimo vettore del ciclo precedente
Ox (1,1) = 0;                                   %Riporta a zero la massa bruciata
Odxb = zeros(siz,1);
for i = 1:siz-1                                  %Secondo ciclo
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
for i = 1:siz-1                                  %Terzo ciclo
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
    %Richiama la legge del manovellismo dalla libreria per conoscere V e dv
    %al variare di teta
[V(i),dv(i)] = OdLibv2.Vdv (fgr,teta,vvc,vvs,rsul);
%Usa la massa, la temperatura e il volume per calcolare la pressione
%tramite la legge dei gas ideali
 P (i)= (r*Ox(i,3)*Ox(i,2)/V(i))/100000;

end

    % Ciclo di pompaggio
    
    Pdue = [P P];
    Vdue = [V V];
    
    
% Valuto il lavoro utile indicato 

Lui = zeros (siz,1);                            %Inizializzo 

for i = 1:siz-1
    
Lui(i+1) = Lui(i) + P(i)*100000 * dv(i) * dteta;

end



Lui = Lui/1000;                                     %Lavoro utile indicato

% Calcolo i dati tecnici del motore
    
Lue = Lui * etam;                                   %Lavoro effettivo                   
display 'Potenza in Kw'
potkw = nc * Lue(siz) * rpm / 60 / 2                %Potenza in Kw
display 'Coppia'
copnm = potkw * 1000 * 60 / 2 / pi / rpm            %Coppia     
display 'Pressione media'
pmi = Lui(siz) / vvs /100                           %Pressione media
mtrap = Ox (siz,3);
mf = mtrap / (alfa + 1);
display 'Consumo specifico'
cs = (mf * 1000) / (Lue(siz) / 3600) / etab         %Consumo specifico 
display 'Riempimento'
lambdav = mtrap / mref                              %Riempimento 



 Odxb= Odxb*100000/4;                                            %Aggiustamento di scala 
  
  % Salva i dati su file esterno

%   save Calcolo.mat;
   
  % plot
  
  figure (5)                                                     %Legge di rilascio del calore
  plot (tetaciclo(250:1000),Ox(250:1000,1)/Ox(250:1000,3)); 
 hold on                                                        %Velocità di rilascio
  plot (tetaciclo(250:1000),Odxb(250:1000,1));
  
 figure (7)                                             %Pano P-V
  plot(V,P); 
  
   figure (8)                                                   %Andamento del lavoro indicato durante il ciclo
  plot (tetaciclo,Lui); 
  
   figure (2)                                 %Ciclo P-teta
  plot (tetaciclo,P);  
%     
%     subplot (3,1,3)
  figure (3)                                %Legge di variazione del volume 
  
  
    plot (tetaciclo,V)
    
    % Ciclo di pompaggio
    
    Pdue = [P P];
    Vdue = [V V];
    
    figure (10)
    
    plot ( Vdue (1700:3800),Pdue(1700:3800)  )
    
     figure (4)                                         %Grafico T-teta
  plot (tetaciclo,Ox(:,2));
  
 figure (6)                                         %Grafico massa evolvente del cilindro 
  plot (tetaciclo,Ox(:,3));
  

    
