% MOFlowEarthCMB
%% NEED TO CHECK PE CALCS

% Using the Th-U decay chain
% lines that need to be changed for new planets are marked with ***
% **** change coretemp in cool2clement to make smooth profile

% code that integrates atmospheric growth and time to solidification with solid phase calculations
% calculates F(j) from planet using T_surf(j-1) and emissivity(i-1)
% then pressure of atmosphere (j)
% then time since last step (j)
% then then saturation levels of CO2 and H2O in magma at surface (j)
% degasses interior for delta(j) time using velocity scaling law (j)
% then calculates emissivity (j)
% adiabatically continues interior temperature to bottom of viscous boundary: this is T_top(j)
% uses F(j) to calculate delta(T) across boundary and thus obtain new
% T_surf(j)

close all;
%clear all;

H2Oliquid(1) = 0.00; %0.5;    %0;  % in mass percent
CO2liquid(1) =  0.00;%0.1;%0.6; %

DM_string = num2str(DM/1000);

% initial values for calculations in this program
CMB = 3470000;                  % *** m, radius of core-mantle boundary
R = 6378000;                    % *** m, total radius of planet
%RM = (R - CMB); %2000000;                   % *** m, depth of magma ocean
g = 9.8;                        % *** m/sec2
adiabslope = 0.33/1000;         %*** K/m, slope of adiabat
tfinal = 50*3.14e13;            % *** sec total time of conductive cooling in cool2clement
tempcore = 1600 + 273;         % *** in K (2100C = 2373 K; 1900C = 2173K) match to ending T in MOFlow
InitP = -0.0374*((R-DM)/1000) + 238.5372;   % ***[GPa] at bottom of MO
Tsolidend = 600;                %*** temp when all interior is solid according to our calcs; used to distribute latent heat
Tsurflatent = 1500;             %*** surface temp below which latent heat begins to be phased out
Mantlemass = (R^3 - (R - DM)^3)/(R^3 - CMB^3)*4.032e+024; %*** kg, mass of MO
%sum
name = ([num2str(H2Oliquid(1)),'% H_2O, ',num2str(CO2liquid(1)),'% CO_2']);

maxstep = 998;  % number of steps = number of volume fractions
tage = 4.56*10^9*3.14*10^7; % age of planet, in sec
kth = 3;            % W/mK
H = 418700;         % J/kg
Cp = 1256.1;        % J/kg K
cpcore = 800;               % J/kg K
sigma = 5.67e-8;    % J/m2K4sec (Boltzmann's constant)
kappa = 1e-6;       % m2/sec (can calculate using Zahnle et al 1988 Ra formulation)
rho = 3000;         % kg/m3
rhocore = 7500;             % kg/m^3
solidrho = 4000;    % kg/m3
alpha = 3e-5;       % K-1
r = zeros(1, maxstep);
r(1) = (R-DM);
eta = 1;            % Pas
kwater = 0.01;      % m2/kg Yamamoto '52 for water for emissivity calculations
kcarbon = 0.01; %0.001; %0.05;     % m2/kg Pujol and North 2003
mwater = 18e-3;     % kg/mol
mcarbon = 44e-3;    % kg/mol
po = 101325;        % Pa
Mantlevolume = (4/3)*pi*((R)^3 - (r(1))^3); % m3 volume of mantle
%Tsolid(1) = (-1.1601e-007)*(r(1)/1000)^3 + 0.0014*(r(1)/1000)^2 +
%-6.3821*(r(1)/1000) + 1.4439e+004;  - old
%Tsolid(1) = -0.000000000086301e3*(r(1)/1000)^3 + 0.000000887643612e3*(r(1)/1000)^2 +...
%    (-0.003265497414172e3*(r(1)/1000)) + 8.310438900359062e3;

solidus = makeSolidus();

Tsolid(1) = solidus(RtoP(r(1)));

Tsurf(1) = Tsolid(1) - adiabslope*(R - r(1)) - H/Cp;        % continue up along adiabat
liqfrac(1) = 1;
time(1) = 0; deltatime(1) = 0; veloc(1) = 0;
Hmassflux(1) = 0; Cmassflux(1) = 0; Hatmadd(1) = 0; HMOfactor(1) = 0; Catmadd(1) = 0; CMOfactor(1) = 0;
delt(1) = 0; Cfractionout(1) = 0; Hfractionout(1) = 0;

if H2Oliquid(1) > 0; 
    initialatmospherewater;  
else HPatm(1) = 0; 
    H2Oliquid2(1) = 0; 
end
% this program partitions the initial volatile quantities into atm and MO

if CO2liquid(1) > 0; 
    initialatmospherecarbon;  
else CPatm(1) = 0; 
    CO2liquid2(1) = 0; 
end
% this program partitions the initial volatile quantities into atm and MO

% 1: SiO2, 2: Al2O3, 3: FeO, 4: MgO, 5: CaO, 6: Sm, 7: Nd, 8: Th, 9: U, 10: OH, 11: C
% Earth: Hart and Zindler - Na2O plus 2xchondritic traces plus chondritic water, all mass percent

liq_comp = [45.96 4.06 7.54 37.78 3.21 0.00001472 0.00004524 0.00000294 0.00000081 0 0];
%0.00000243,0.0000104,... % changing from Lu-Hf to Th-U
% U & Th from Anders and Grevesse CI chondrite table, column G.,

% Molar weights of oxides and elements.  Putting this here to make it easy
% to change in once place, rather than each individual mineral file,
% especially if we resort to using Lu/Hf to confirm Carlson's suggestions.
MW = [  60.09,...   % SiO2   (1)   g/mol for oxides
       101.96,...   % Al2O3  (2)
        71.846,...  % FeO    (3)
        40.311,...  % MgO    (4)
        56.077,...  % CaO    (5)
       150.36,...   % Sm     (6) (these in molar ppm after division)
       144.24,...   % Nd     (7)
       232.0381,... % Th     (8)
       238.0289,... % U      (9)
       17.00,...    % OH     (10)   g/mol again
       28.00];      % CO     (11)
 
       %174.967,...  % Lu     (8)
       %178.49,...   % Hf     (9)  
    
liq_comp = 100*liq_comp./(100-(H2Oliquid2(1)+CO2liquid2(1))); % normalized MO oxides so volatiles remain at values specified above
liq_comp(10) = H2Oliquid2(1); liq_comp(11) = CO2liquid2(1);

% initial values for calculations in fractionate.m

% post-perovskite layer above the CMB, as magnesiowustite (sp?) is unstable
% at the relevant temperatures and pressures
intliqall = 0.00; 

intliq0 = intliqall;	% post-perovskite
ppv0 = (1 - intliq0); % layer 0
Layer0P = 120; %

intliq1 = intliqall;	% perovskite and magnesiowustite
per1 = (1 - intliq1)*0.95; % layer 1
mw1 = (1 - intliq1)*0.05;
Layer1P = 22;

intliq2 = intliqall;	% gamma olivine and majorite
gamma2 = (1 - intliq2)*0.45; %0.55; % layer 2
maj2 = (1 - intliq2)*0.55; %0.45;
Layer2P = 18.00;

intliq3 = intliqall;	% beta olivine, majorite, cpx
beta3 = (1 - intliq3)*0.55; %0.5; %0.4; % layer 3
maj3 = (1 - intliq3)*0.2; %0.25; %0.35;
cpx3 = (1 - intliq3)*0.25; %0.25;
Layer3P = 15.00;

intliq4 = intliqall;	% garnet, cpx, opx, olivine
gar4 = (1 - intliq4)*0.10; %0.10; %0.15; % layer 4 garnet settling is removed
cpx4 = (1 - intliq4)*0.20;
opx4 = (1 - intliq4)*0.20; %0.15;
alpha4 = (1 - intliq4)*0.50;
Layer4P = 2.5;

intliq5 = intliqall;	% spinel, opx, cpx, olivine
spin5 = (1 - intliq5)*0.05; %0.10; %0.15;  % layer 5
cpx5 = (1 - intliq5)*0.2; %0.25;
opx5 = (1 - intliq5)*0.20;
alpha5 = (1 - intliq5)*0.55; %0.50; %0.45; %0.40;
Layer5P = 1.0;

intliq6 = intliqall;	% plagioclase, cpx, opx, olivine
plag6 = (1 - intliq6)*0.05; %0.10; %0.15; % layer 6
cpx6 = (1 - intliq6)*0.25;
opx6 = (1 - intliq6)*0.25; %0.20;
alpha6 = (1 - intliq6)*0.45; %0.50; %0.45; %0.40;
Layer6P = 0.6;

intliq7 = intliqall;   % kind of a proxy for all the remaining evolved liquid
plag7 = (1 - intliq7)*0.3; %  no liquid evolution
cpx7 = (1 - intliq7)*0;
opx7 = (1 - intliq7)*0.3;
alpha7 = (1 - intliq7)*0.3;
Layer7P = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P(1) = InitP;           % P at bottom of magma ocean
liquid(1,:) = liq_comp; intliq(1,:) = liq_comp;
MG(1,1) = (liq_comp(4)/40.311)./(liq_comp(4)/40.311 + liq_comp(3)/71.846); % Mg1 is liquid, 2 is solid

% parameters for sortandinvert.m
dFdT = 0.3;                 % weight % melted per degree rise above solidus
Planetmass = 5.9742e24;     %*** mass of Earth, from NASA [kg]
Coremass = 1.9416e024;      %*** mass of core [kg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%begin
Patm(1) = CPatm(1) + HPatm(1); % calculated above when initialatmosphere routines are called
HMatm(1) = (HPatm(1)*4*pi*R^2)/g; CMatm(1) = (CPatm(1)*4*pi*R^2)/g;
Matm(1) = HMatm(1) + CMatm(1);

taustarw = zeros(1,maxstep);
taustarw(1) = ((3*HMatm(1))/(8*pi*R^2))*(((kwater*g)/(3*po))^0.5);  % Abe and Matsui 1986

taustarc = zeros(1,maxstep);
taustarc(1) = ((3*CMatm(1))/(8*pi*R^2))*(((kcarbon*g)/(3*po))^0.5);

emiss = zeros(1,maxstep);
emiss(1) = (2 /(taustarw(1)+taustarc(1) + 2));

vol = zeros(1,maxstep);
vol(1) = 1;

flux = zeros(1,maxstep);

k = 1;  m = 1; % loops for abbreviated data for output - see bottom of file

magnesio_thermal = zeros(maxstep, 2);

%% Vectors for later graphing tests
mass_solidified = 0;

all_liquid_composition = zeros(maxstep, length(liq_comp));
all_liquid_composition(1, :) = liq_comp;

mantle_mass_vector = zeros(maxstep, 1);
mantle_mass_vector(1, :) = mass_solidified;

mantle_mass_by_layer = zeros(maxstep, 1);
mantle_mass_by_layer(1, :) = mass_solidified;

% mantle_volume_by_layer = ones(maxstep, 1);
% mantle_volume_by_layer = mantle_volume_by_layer.*(Mantlevolume./1000);

solid_comp_by_layer = zeros(maxstep, length(liq_comp));
solid_comp_by_layer(1, :) = liq_comp;

residual_liquids_vector = zeros(maxstep, 1);
residual_liquids_vector(1, :) = Mantlemass;

samarium = zeros(maxstep, 1);
samarium_percent = zeros(maxstep, 1);
samarium(1, :) = liq_comp(6)./150.36;
samarium_percent = samarium(1, :)*mass_solidified;

for j = 2:1:maxstep    % each step is one-tenth of a percent solidification by volume
%     display(j)
    vol(j) = j;

    % calculate emissivity, time, saturation levels in magma for i
    taustarw(j) = ((3*HMatm(j-1))/(8*pi*R^2))*(((kwater*g)/(3*po))^0.5);  % Abe and Matsui 1986
    taustarc(j) = ((3*CMatm(j-1))/(8*pi*R^2))*(((kcarbon*g)/(3*po))^0.5);
    emiss(j) = (2 /(1.5*taustarw(j)+0.5*taustarc(j) + 2));    % 1.5 and 0.5 come from relative absorption wavelength widths

    % calculate F
    flux(j) = emiss(j)*sigma*((Tsurf(j-1)+273)^4 - 250^4);

    liqfrac(j) = (1000 - j)/1000;    % depends on interpretation of j as 0.25% solidification
    delr(j) = ((R)^3 - (3/(4*pi)*Mantlevolume*(liqfrac(j))))^(1/3) - (r(j-1)); %m, thickness of new solidified layer; depends on interpretation of j as 0.25% solidification
    r(j) = r(j-1) + delr(j);  % r [m]

    % continue T to bottom of viscous boundary layer to get Ttop (T's following solidification of this step)
    %delt(j) = 10*[(rho*alpha*g*flux(j))/(eta*kappa*k)]^(-1/4); XX what is k here?
    %Tsolid(j) = (-1.1601e-007)*(r(j)/1000)^3 + 0.0014*(r(j)/1000)^2 + -6.3821*(r(j)/1000) + 1.4439e+004; % doesn't seem to need this + (- 10*(0.2*(liqfrac(j))+0.02)^(-1));
%     Tsolid(j) =  -0.000000000086301e3*(r(j)/1000)^3 + 0.000000887643612e3*(r(j)/1000)^2 +...
%         (-0.003265497414172e3*(r(j)/1000)) + 8.310438900359062e3;
    Tsolid(j) = solidus(RtoP(r(j)));
    Ttop(j) = Tsolid(j) - adiabslope*(R - r(j)) - H/Cp;        % continue up along adiabat                                                                                                                                    % depends on interpretation of j as 0.25% solidification
    Tsurf(j) = Ttop(j); % + [flux(j)/((0.1*kthXX??*rho*alpha*g)/(eta*kappa))]^(3/4);      % converted j (volume) to radius to fit equation
    if Tsurf(j-1) < Tsurflatent
        Tsurf(j) = Tsolid(j) - adiabslope*(R - r(j)) ...
            - [(Tsurf(j-1)-Tsolidend)/(Tsurflatent - Tsolidend)]*H/Cp;        % continue up along adiabat                                                                                                                                    % depends on interpretation of j as 0.25% solidification
    end

    % time
    heatterm = (H/Cp)*(r(j)/R)^2 + (R/(3*delr(j)))*(1 - (r(j)/R)^3)*(Tsolid(j-1) - Tsolid(j));
    timerate(j) = flux(j)/(heatterm*rho*Cp);    % solidification rate in m/sec
    deltatime(j) = (delr(j))/timerate(j);   % time to do this increment in sec
    time(j) = time(j-1) + deltatime(j);

    % degas magma ocean
    veloc(j) = (0.6/(0.6)^(1/3))*[(alpha*g*flux(j)*0.5)/(rho*Cp)]^(1/3);
    Hmassflux(j) = (liqfrac(j)*Mantlemass);
    Cmassflux(j) = (liqfrac(j)*Mantlemass);

    % calculate saturation values in magma
    H2Osaturation(j) = 100*[(2.08e-6)*(HPatm(j-1)).^0.52] + 0.3;  % H2O saturation in mass percent, P in Pa
    CO2saturation(j) = 100*[(2.08e-6)*(CPatm(j-1))^0.45] + 0.05;  % CO2 saturation in mass percent, P in Pa

    % call program to fractionate solid phases and send back the amount of volatiles to remove from the liquid
    fractionatedeepEarth

    %%% water degassing routines
    if liquid(j,10) < H2Osaturation(j); Hatmadd(j) = 0; HMOfactor(j) = 0;
    else
        Hatmadd(j) = Hmassflux(j)*((liquid(j,10))- H2Osaturation(j))/100;  % kg of water into atmosphere
        HMOfactor(j) = Hmassflux(j)/(liqfrac(j)*Mantlemass);          % fraction of MO that is degassed in this step
        liquid(j,10) = (HMOfactor(j)*H2Osaturation(j) + (1-HMOfactor(j))*liquid(j-1,10));
    end

    %%% carbon degassing routines
    if liquid(j,11) < CO2saturation(j); Catmadd(j) = 0; CMOfactor(j) = 0;
    else
        Catmadd(j) = Cmassflux(j)*((liquid(j,11))- CO2saturation(j))/100;    % kg of carbon into atmosphere
        CMOfactor(j) = Cmassflux(j)/(liqfrac(j)*Mantlemass);            % fraction of MO that is degassed in this step
        liquid(j,11) = (CMOfactor(j)*CO2saturation(j) + (1-CMOfactor(j))*liquid(j-1,11));
    end

    liquid(j,:) = 100*liquid(j,:)./(sum(liquid(j,:)));  % renormalized liquid composition

    % calculate new atmospheric mass and pressure
    HMatm(j) = HMatm(j-1) + Hatmadd(j);
    CMatm(j) = CMatm(j-1) + Catmadd(j);
    Matm(j) = Matm(j-1) + Hatmadd(j) + Catmadd(j);

    Patm(j) = ((CMatm(j-1) + HMatm(j-1))*g)/(4*pi*R^2);  %GPa
    HPatm(j) = ((HMatm(j-1))*g)/(4*pi*R^2);   % translates P in GPa to Matm in kg
    CPatm(j) = ((CMatm(j-1))*g)/(4*pi*R^2);

    % calculate fraction of initial volatile degassed; if loop prevents division by zero
    if CO2liquid(1)>0; Cfractionout(j) = (CMatm(j))/(Mantlemass*CO2liquid(1)/100); else Cfractionout(j)=0;end
    if H2Oliquid(1)>0; Hfractionout(j) = (HMatm(j))/(Mantlemass*H2Oliquid(1)/100); else Hfractionout(j)=0;end

    % fill in a few final numbers
    MG(1,2) = MG(2,2);
    solid(1,:) = solid(2,:);
    D(1) = D(2);
    Do(1) = Do(2);
    Dsol(1) = Dsol(2);

%     % for a figure
%     if rem(j,100) == 0
%         tfig(k) = time(j)/3.1536e13;
%         Tsurffig(k) = Tsurf(j);
%         Tsolidfig(k) = Tsolid(j);
%         H2Ofig(k) = liquid(j,10);
%         CO2fig(k) = liquid(j,11);
%         rfig(k) = r(j);
%         k = k+1;
%     end
% 
%     if j == 2
%         tfig(k) = time(j)/3.1536e13;
%         Tsurffig(k) = Tsurf(j);
%         Tsolidfig(k) = Tsolid(j);
%         H2Ofig(k) = liquid(j,10);
%         CO2fig(k) = liquid(j,11);
%         rfig(k) = r(j);
%         k = k+1;
%     end

end

%figure(1);title(['Density with depth for model:  ', name]); hold on; plot(D,r/1000, 'r'); xlabel('density [kg/m3]'); ylabel('radius, km');
%figure(2);title(['Reference density with depth for model:  ', name]); hold on; plot(Do,r/1000, 'r'); xlabel('density at 1 atm and 1 deg C [kg/m3]'); ylabel('radius, km');
%figure(3); %title(['Reference density with depth for model:  ', name]); 
 %   hold on; 
    %plot(Dsol,r/1000, 'r'); 
    %xlabel('density at 1 atm and solidus temperature [kg/m^3]'); 
    %ylabel('radius, km');

sortandinvertperovskite
%%%calculate new emissivity with volatile release from melting
taustarwend = ((3*(HMatm(j)+Hnewatm))/(8*pi*R^2))*(((kwater*g)/(3*po))^0.5);  % Abe and Matsui 1986
taustarcend = ((3*(CMatm(j)+Cnewatm))/(8*pi*R^2))*(((kcarbon*g)/(3*po))^0.5);
endemiss = (2 /(taustarwend+taustarcend + 2));    % H2O/CO2 relative absorption wavelength widths =1.5/0.5 not considered here
%%%%%%%%%%%%
%cool2clementwhole
graphsdeep
%leanGraphs
Ddoubleprime

epNd