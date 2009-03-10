function density = majoritedensity(Mgnum, Perc_Ca, Perc_AlMgFe, Perc_MgFe, P, T)
 
% calculates density of majorite
% calculate Vo from Vo298 and integral for alpha dT
% calculate V(P,T) from Birch-Murnaghan equation at given P (molar volumes at P and T)
% using KT and KpT, each calculated from their own derivatives
% calculate density from V(P,T) and molar weights of each variety of majorite
% physical parameters from Akaogi et al (2002), Hama and Suita (2001), Bertka and Fei (1998)
 
% at first, physical parameter matrix is size (7 x 3)
% where each row is a parameter
% in order 1 = Vo298, 2 = KT298, 3 = dK/dP, 4 = dK/dT, 5 = ao, 6 - a1, 7 = a2
% and each column is a phase: 1 = CaAl, 2 = MgFeAl, 3 = MgFe
 
% below parameters are added: 8 = Vo(T), 9 = KT, 10 = V(P,T) in m3/mole
%                               11 = kg/mol of formula units of phase, 12 = density of phase in kg/m3
%                               13 = calculated alpha, /deg
 
% physical parameters for Ca majorite (Ca)3(Al)2(Si)3(O)12
par(1,1) = 125.12;      % Vo298, cm3/mole
par(2,1) = 168;         % KT298, GPa
par(3,1) = 6.2;         % dK/dP
par(4,1) = -0.022;      % dK/dT
par(5,1) = 1.951e-5;    % ao
par(6,1) = 8.089e-9;    % a1
par(7,1) = -4.972e-1;   % a2
par(11,1) = (3*56.077  + 101.96 + 3*60.09)/1000; % kg/moles of formula
 
% physical parameters for Al, Mg, Fe majorite (Mg, Fe)3(Al)2(Si)3(O)12
par(1,2) = Mgnum*113.08 + (1-Mgnum)*115.43;         % Vo298, cm3/mole
par(2,2) = Mgnum*179 + (1-Mgnum)*175; %261;               % KT298, GPa
par(3,2) = 4.00;                                    % dK/dP
par(4,2) = -0.022;                                  % dK/dT
par(5,2) = Mgnum*2.311e-5 + (1-Mgnum)*1.776e-5;     % ao
par(6,2) = Mgnum*5.956e-9 + (1-Mgnum)*1.214e-8;     % a1
par(7,2) = Mgnum*(-0.4538) + (1-Mgnum)*(-0.5071);   % a2
par(11,2) = (3*(Mgnum*40.311 + (1-Mgnum)*71.846) + 101.96 + 3*60.09)/1000;  % kg/moles of formula
 
% physical parameters for Mg, Fe majorite (Mg,Fe)4(Si)4(O)12
par(1,3) = 113.99; %Mgnum*85.5 + (1-Mgnum)*88.2; % Vo298, cm3/mole
par(2,3) = 160; %161;             % KT298, GPa
par(3,3) = 4.0;             % dK/dP
par(4,3) = -0.022;          % dK/dT
par(5,3) = 2.874e-5;        % ao
par(6,3) = 2.886e-9;        % a1
par(7,3) = -0.5443;         % a2
par(11,3) = (4*(Mgnum*40.311 + (1-Mgnum)*71.846) + 4*60.09)/1000;  % kg/moles of formula
 
for i = 1:3     % loop through each phase
 
      par(13,i) = (par(5,i)*(T+298) + (par(6,i)/2)*((T+298)^2) - par(7,i)*(1/(T+298)))...
         - (par(5,i)*(298) + (par(6,i)/2)*(298^2) - par(7,i)*(1/298));
      par(8,i) = par(1,i)*exp(par(13,i));        % Vo(T)
      par(9,i) = par(2,i) + par(4,i)*(T+298);    % KT(T)
      Kot = par(9,i);                            % for consistency in function BirchMurnsolid
      Kpt = par(3,i);
      Vo = par(8,i);
      VP(i) = BirchMurnsolid(Kot, Kpt, Vo, P);
      par(10,i) = VP(i)/1e6;                     % convert molar volume of m3/mol
      par(12,i) = par(11,i)/par(10,i);           % phase density in kg/m3
 
end
 
density = Perc_Ca*par(12,1) + Perc_AlMgFe*par(12,2) + Perc_MgFe*par(12,3);
