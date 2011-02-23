function density = clinopyroxenedensity(Mgnum, CaMg, P, T)
 
% calculates density of clinopyroxene
% calculate Vo from Vo298 and integral for alpha dT
% calculate V(P,T) from Birch-Murnaghan equation at given P (molar volumes at P and T)
% using KT and KpT, each calculated from their own derivatives
% calculate density from V(P,T) and molar weights of each variety of pyroxene
% physical parameters from Bertka and Fei (1998) 
 
% at first, physical parameter matrix is size (7 x 3)
% where each row is a parameter
% in order 1 = Vo298, 2 = KT298, 3 = dK/dP, 4 = dK/dT, 5 = ao, 6 - a1, 7 = a2 
% below parameters are added: 8 = Vo(T), 9 = KT, 10 = V(P,T) in m3/mole
%                               11 = kg/mol of formula units of phase, 12 = density of phase in kg/m3
%                               13 = calculated alpha, /deg

Mg = Mgnum/(Mgnum + (1-Mgnum) +CaMg);
Fe = (1-Mgnum)/(Mgnum + (1-Mgnum) +CaMg);
Ca = CaMg/(Mgnum + (1-Mgnum) +CaMg);

% physical parameters for Ca pyroxene (Mg,Fe,Ca)2(Si)2(O)6

par(1,1) = Mgnum*(CaMg*66.04 + (1-CaMg)*62.70) + (1-Mgnum)*(CaMg*67.87 + (1-CaMg)*65.89); %Mgnum*66.04 + (1-Mgnum)*67.87;       % Vo298, cm3/mole 
par(2,1) = Mgnum*(CaMg*113 + (1-CaMg)*107) + (1-Mgnum)*Mgnum*(CaMg*119 + (1-CaMg)*101);   % KT298, GPa
par(3,1) = Mgnum*(CaMg*4.8 + (1-CaMg)*10.2) + (1-Mgnum)*(CaMg*4.2 + (1-CaMg)*4.2);        % dK/dP
par(4,1) = Mgnum*(CaMg*-0.02 + (1-CaMg)*-0.00015) + (1-Mgnum)*(CaMg*-0.02 + (1-CaMg)*0);   % dK/dT
par(5,1) = Mgnum*(CaMg*3.33e-5 + (1-CaMg)*5.05e-5) + (1-Mgnum)*(CaMg*2.98e-5 + (1-CaMg)*6.32e-5);   % ao
par(6,1) = 2.694e-9;               % a1
par(7,1) = -0.5588;               % a2
par(11,1) = (2*(Mg*40.311 + Fe*71.846 + Ca*56.077) + 2*60.09)/1000;  % kg/moles of formula


par(13,1) = (par(5,1)*(T+298) + (par(6,1)/2)*((T+298)^2) - par(7,1)*(1/(T+298)))...
        - (par(5,1)*(298) + (par(6,1)/2)*(298^2) - par(7,1)*(1/298));
par(8,1) = par(1,1)*exp(par(13,1));        % Vo(T)
par(9,1) = par(2,1) + par(4,1)*(T+298);    % KT(T)
Kot = par(9,1);                            % for consistency in function BirchMurnsolid
Kpt = par(3,1);
Vo = par(8,1);
%fprintf('Mgnum=%.2f, CaMg=%.2f, P=%.2f, T=%.2f \n', Mgnum, CaMg, P, T)
% if Mgnum < 0.29 
%     fprintf('clino Mg num < 0.29\n')
if Mgnum < 0.01
    fprintf('clino Mg num < 0\n')
elseif Mgnum > 1
    fprintf('clino Mg num > 1\n')
end
VP(1) = BirchMurnsolid(Kot, Kpt, Vo, P);
par(10,1) = VP(1)/1e6;                     % convert molar volume of m3/mol
par(12,1) = par(11,1)/par(10,1);           % phase density in kg/m3

density = par(12,1);
%fprintf('clino density =%.2f \n', density)
