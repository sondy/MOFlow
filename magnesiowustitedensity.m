function density = magnesiowustitedensity(Mgnum, P, T)

% calculates density of magnesiowustite
% calculate Vo from Vo298 and integral for alpha dT
% calculate V(P,T) from Birch-Murnaghan equation at given P (molar volumes at P and T)
% using KT and KpT, each calculated from their own derivatives
% calculate density from V(P,T) and molar weights of each variety of perovskite
% physical parameters from Fabrichnaya (1999); Akaogi et al (2002), Hama and Suita (2001), Bertka and Fei (1998)

% at first, physical parameter matrix is size (7 x 2)
% where each row is a parameter
% in order 1 = Vo298, 2 = KT298, 3 = dK/dP, 4 = dK/dT, 5 = ao, 6 - a1, 7 = a2
% and each column is a phase: 1 = Fe, 2 = Mg

% below parameters are added: 8 = Vo(T), 9 = KT, 10 = V(P,T) in m3/mole
%                               11 = kg/mol of formula units of phase, 12 = density of phase in kg/m3
%                               13 = calculated alpha, /deg; 14 = a3 (for four-term Fabrichnaya calculation)

% physical parameters for FeO magnesiowustite
Vol298(1) = 12.25; %18.4; %12.25;       % Vo298, cm3/mole
KT298(1) = 158;         % KT298, GPa
Kpt(1) = 4.13;        % dK/dP
par(4,1) = -0.031;           % dK/dT
a0(1) = 0.355e-4;    % ao
a1(1) = 1.1e-8;    % a1
a2(1) = 0;   % a2
rho_mol(1) = 71.846/1000;  % kg/moles of formula

% physical parameters for MgO magnesiowustite
Vol298(2) = 11.25; %18.4; %11.25; % Vo298, cm3/mole
KT298(2) = 160; % KT298, GPa
Kpt(2) = 4.13; % dK/dP
par(4,2) = -0.0272; % dK/dT
a0(2) = 0.3768e-4; % ao
a1(2) = 0.7404e-8;    % a1
a2(2) = -0.7445; % a2
rho_mol(2) = 40.311/1000;  % kg/moles of formula

VP = zeros(1, 2);

for i = 1:2     % loop through each phase
%    display(j);
%    display(i);
    par(13,i) = (a0(i)*(T+298) + (a1(i)/2)*((T+298)^2) - a2(i)*(1/(T+298)))...
        - (a0(i)*(298) + (a1(i)/2)*(298^2) - a2(i)*(1/298));
    Vo(i) = Vol298(i)*exp(par(13,i));            % Vo(T)
    Kot(i) = KT298(i) + par(4,i)*(T+298);        % KT(T)
    
%     magnesio_thermal(j, i) = Kot( i);
%     Kot = Kot(i);                                 % for consistency in function BirchMurnsolid
%     Kpt = Kpt(i);
    Vo = Vo(i);
    VP(i) = BirchMurnsolid(Kot, Kpt, Vo, P);
    molvol(i) = VP(i)/1e6;              % convert molar volume of m3/mol
    phase_rho(i) = rho_mol(i)/molvol(i);    % phase density in kg/m3

end

density = (1-Mgnum)*phase_rho(1) + (Mgnum)*phase_rho(2);

