function density = spineldensity(Mgnum, P, T)
 
% calculates density of spinel
% calculate Vo from Vo298 and integral for alpha dT
% calculate V(P,T) from Birch-Murnaghan equation at given P (molar volumes at P and T)
% using KT and KpT, each calculated from their own derivatives
% calculate density from V(P,T) and molar weights of each variety of garnet
% physical parameters from Akaogi et al (2002), Hama and Suita (2001), Bertka and Fei (1998) 
 
% at first, physical parameter matrix is size (7 x 3)
% where each row is a parameter
% in order 1 = Vo298, 2 = KT298, 3 = dK/dP, 4 = dK/dT, 5 = ao, 6 - a1, 7 = a2
% and each column is a phase: 1 = CaAl, 2 = MgFeAl, 3 = MgFe
 
% below parameters are added: 8 = Vo(T), 9 = KT, 10 = V(P,T) in m3/mole
%                               11 = kg/mol of formula units of phase, 12 = density of phase in kg/m3
%                               13 = calculated alpha, /deg
 
% physical parameters for spinel (Fe, Mg)Al2O4
par(1,1) = Mgnum*39.78 + (1-Mgnum)*40.75;       % Vo298, cm3/mole
par(2,1) = Mgnum*194.5 + (1-Mgnum)*212.0;         % KT298, GPa
par(3,1) = 4.00;        % dK/dP
par(4,1) = -0.0001;      % dK/dT
par(5,1) = Mgnum*4.31e-5 + (1-Mgnum)*3.95e-5;    % ao
par(6,1) = 0; %8.089e-9;    % a1
par(7,1) = 0; %-4.972e-1;   % a2
par(11,1) = ((Mgnum*40.311 + (1-Mgnum)*71.846) + 101.96)/1000; % kg/moles of formula

     par(13,1) = (par(5,1)*(T+298) + (par(6,1)/2)*((T+298)^2) - par(7,1)*(1/(T+298)))...
        - (par(5,1)*(298) + (par(6,1)/2)*(298^2) - par(7,1)*(1/298));
     par(8,1) = par(1,1)*exp(par(13,1));        % Vo(T)
     par(9,1) = par(2,1) + par(4,1)*(T+298);    % KT(T)
     Kot = par(9,1);                            % for consistency in function BirchMurnsolid
     Kpt = par(3,1);
     Vo = par(8,1);
     VP(1) = BirchMurnsolid(Kot, Kpt, Vo, P);
     par(10,1) = VP(1)/1e6;                     % convert molar volume of m3/mol
     par(12,1) = par(11,1)/par(10,1);           % phase density in kg/m3

density = par(12,1);
