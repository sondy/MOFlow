function density = postperovskitedensity(Mgnum, CaMg, P, T)

CaMg=0; % Not set up for Ca in ppv
% calculates density of postperovskite (Mg,Fe)SiO3
% calculate Vo from Vo298 and integral for alpha dT
% calculate V(P,T) from Birch-Murnaghan equation at given P (molar volumes at P and T)
% using KT and KpT, each calculated from their own derivatives
% calculate density from V(P,T) and molar weights of each variety of pyroxene
% physical parameters from Akaogi et al (2002), Hama and Suita (2001), Bertka and Fei (1998) 
 
% at first, physical parameter matrix is size (7 x 3)
% where each row is a parameter
% in order 1 = Vo298, 2 = KT298, 3 = dK/dP, 4 = dK/dT, 5 = ao, 6 - a1, 7 = a2 
% below parameters are added: 8 = Vo(T), 9 = KT, 10 = V(P,T) in m3/mole
%                               11 = kg/mol of formula units of phase, 12 = density of phase in kg/m3
%                               13 = calculated alpha, /deg

Mg = Mgnum/(Mgnum + (1-Mgnum) +CaMg);   
Fe = (1-Mgnum)/(Mgnum + (1-Mgnum) +CaMg);

% physical parameters for postperovskite (Mg,Fe)SiO3

par(1,1) = 24.66; % Vo298, cm3/mole 
par(2,1) = 219;   % KT298, GPa
par(3,1) = 4.0;        % dK/dP
par(4,1) = -0.015;   % dK/dT
par(5,1) = 1.48e-5;   % ainf (Guignot 2007)
par(6,1) = 396;               % b
par(7,1) = 0;               % not used
par(11,1) = ((Mg*40.311 + Fe*71.846) + 60.09)/1000;  % kg/moles of formula

   par(13,1) = par(5,1).*exp(-par(6,1)/(1+(T+298)));
   par(8,1) = par(1,1)*exp(par(13,1));        % Vo(T)
   par(9,1) = par(2,1) + par(4,1)*(T+298);    % KT(T)
   Kot = par(9,1);                            % for consistency in function BirchMurnsolid
   if Kot < 0
       PrintCaller
   end
   Kpt = par(3,1);
   Vo = par(8,1);
   VP(1) = BirchMurnsolid(Kot, Kpt, Vo, P);
   par(10,1) = VP(1)/1e6;                     % convert molar volume of m3/mol
   par(12,1) = par(11,1)/par(10,1);           % phase density in kg/m3

density = par(12,1);
