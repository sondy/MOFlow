function density = olivinedensity(Mgnum, P, T)
 
% calculates density of olivine
% calculate Vo from Vo298 and integral for alpha dT
% calculate V(P,T) from Birch-Murnaghan equation at given P (molar volumes at P and T)
% using KT and KpT, each calculated from their own derivatives
% calculate density from V(P,T) and molar weights of each variety of olivine
% physical parameters from Akaogi et al (2002), Hama and Suita (2001), Bertka and Fei (1998) 
 
% at first, physical parameter matrix is size (7 x 3)
% where each row is a parameter
% in order 1 = Vo298, 2 = KT298, 3 = dK/dP, 4 = dK/dT, 5 = ao, 6 - a1, 7 = a2 
% below parameters are added: 8 = Vo(T), 9 = KT, 10 = V(P,T) in m3/mole
%                               11 = kg/mol of formula units of phase, 12 = density of phase in kg/m3
%                               13 = calculated alpha, /deg
 
% physical parameters for (Mg, Fe) olivine (Mg, Fe)2(Si)1(O)4
par(1,1) = Mgnum*43.6 + (1-Mgnum)*46.29;       % Vo298, cm3/mole 
par(2,1) = Mgnum*129 + (1-Mgnum)*137.9;           % KT298, GPa
par(3,1) = Mgnum*5.37 + (1-Mgnum)*4.0;           % dK/dP
par(4,1) = Mgnum*(-0.0224) + (1-Mgnum)*(-0.0258); % dK/dT
par(5,1) = Mgnum*3.034e-5 + (1-Mgnum)*2.386e-5; % ao
par(6,1) = Mgnum*7.422e-9 + (1-Mgnum)*11.53e-9; % a1
par(7,1) = Mgnum*(-0.5381) + (1-Mgnum)*(-0.0518);   % a2
par(11,1) = (2*(Mgnum*40.311 + (1-Mgnum)*71.846) + 60.09)/1000;  % kg/moles of formula
 
 
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
