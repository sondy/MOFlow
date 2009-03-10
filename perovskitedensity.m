function density = perovskitedensity(Mgnum, Perc_Al, Perc_Ca, Perc_MgFe, P, T)

% calculates density of perovskite
% calculate Vo from Vo298 and integral for alpha dT
% calculate V(P,T) from Birch-Murnaghan equation at given P (molar volumes at P and T)
% using KT and KpT, each calculated from their own derivatives
% calculate density from V(P,T) and molar weights of each variety of perovskite
% physical parameters from Akaogi et al (2002), Hama and Suita (2001), Bertka and Fei (1998)

% does it make sense that Ca=pv is the least dense?


% at first, physical parameter matrix is size (7 x 3)
% where each row is a parameter
% in order 1 = Vo298, 2 = KT298, 3 = dK/dP, 4 = dK/dT, 5 = ao, 6 - a1, 7 = a2
% and each column is a phase: 1 = Al, 2 = Ca, 3 = non-Al, non-Ca

% below parameters are added: 8 = Vo(T), 9 = KT, 10 = V(P,T) in m3/mole
%                               11 = kg/mol of formula units of phase, 12 = density of phase in kg/m3
%                               13 = calculated alpha, /deg

% physical parameters for Al perovskite (Mg,Fe)3(Al)2(Si)3(O)12
par(1,1) = 99.52;       % Vo298, cm3/mole  Akaogi: 99.52; Hama: 160.7
par(2,1) = 261; %234;         % KT298, GPa from Daniel et al. (2001)
par(3,1) = 4.2;         % dK/dP
par(4,1) = 0;           % dK/dT
par(5,1) = 1.982e-5;    % ao
par(6,1) = 8.180e-9;    % a1
par(7,1) = -4.740e-1;   % a2
par(11,1) = (3*(Mgnum*40.311 + (1-Mgnum)*71.846) + 101.96 + 3*60.09)/1000;  % kg/moles of formula

% physical parameters for Ca perovskite (Ca)(Si)(O)3
par(1,2) = 27.45;       % Vo298, cm3/mole  Fei: 27.32*3 = 81.96; Hama: 45*3 = 135
par(2,2) = 236;         % KT298, GPa (286?)
par(3,2) = 3.99;        % dK/dP
par(4,2) = -0.028;      % dK/dT
par(5,2) = 3.2e-5;    % ao
par(6,2) = 9.421e-9;    % a1
par(7,2) = -0.3271;   % a2
par(11,2) = (56.077 + 60.09)/1000;  % kg/moles of formula

% physical parameters for non-Al, non-Ca perovskite (Mg,Fe)(Si)(O)3
par(1,3) = Mgnum*24.5 + (1-Mgnum)*25.6;% Vo298, cm3/mole Fei: Mgnum*73.5 + (1-Mgnum)*76.8;
par(2,3) = Mgnum*261 + (1-Mgnum)*287;         % KT298, GPa
par(3,3) = 4.0;         % dK/dP
par(4,3) = -0.028;      % dK/dT
par(5,3) = 1.982e-5;    % ao
par(6,3) = 8.180e-9;    % a1
par(7,3) = -4.740e-1;   % a2
par(11,3) = ((Mgnum*40.311 + (1-Mgnum)*71.846) + 60.09)/1000;  % kg/moles of formula

VP = zeros(1, 3);

for i = 1:3;     % loop through each phase

    par(13,i) = (par(5,i)*(T+298) + (par(6,i)/2)*((T+298)^2) - par(7,i)*(1/(T+298)))...
        - (par(5,i)*(298) + (par(6,i)/2)*(298^2) - par(7,i)*(1/298));
    %     par(13,i) = par(5,i)*(T+298) - par(5,i)*(298);
    par(8,i) = par(1,i)*exp(par(13,i));        % Vo(T)
    par(9,i) = par(2,i) + par(4,i)*(T+298);    % KT(T)
    Kot = par(9,i);                            % for consistency in function BirchMurnsolid

    if Kot < 1
        PrintCaller
    end

    Kpt = par(3,i);
    Vo = par(8,i);
    VP(i) = BirchMurnsolid(Kot, Kpt, Vo, P);
    par(10,i) = VP(i)/1e6;                     % convert molar volume of m3/mol
    par(12,i) = par(11,i)/par(10,i);           % phase density in kg/m3

end

density = Perc_Al*par(12,1) + Perc_Ca*par(12,2) + Perc_MgFe*par(12,3);
