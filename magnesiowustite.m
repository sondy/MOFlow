function [Mgnum, eqmin] = magnesiowustite(liq_comp)
% MAGNESIOWUSTITE calculates equilibrium magnesiowustite
% (Mg, Fe)O
% liq_comp is in wt%
% MW is the molar weights of oxides
 
FeMgKd = 2;  % Kd = (Fe/Mg)min / (Fe/Mg)liq 2.00 from Fei analysis

M_liq_comp = liq_comp./MW;

load kd_magnesiowustite.dat;
kd_m = kd_magnesiowustite;

KSm = kd_m(1);% 0.001; %0.2;
KNd = kd_m(2);%0.001; %0.06;
KTh = kd_m(3);%0.01; % Change me
KU  = kd_m(4);%0.001; % Change me
KOH = kd_m(5);%0.008;    % KOH = (OH)min/(OH)liq
KCO = kd_m(6);%0.0005;    % complete guesses

% use molar weights and Kds to calculate ratios of elements in mineral
 
FeMgliq = M_liq_comp(3)/M_liq_comp(4);
 
FeMgmin = FeMgKd*FeMgliq;
 
M_eqmin(3) = FeMgmin;
M_eqmin(4) = 1;
M_eqmin(1) = 0;
M_eqmin(2) = 0;
M_eqmin(5) = 0;
M_eqmin(6) = KSm*M_liq_comp(6);
M_eqmin(7) = KNd*M_liq_comp(7);
M_eqmin(8) = KTh*M_liq_comp(8);
M_eqmin(9) = KU*M_liq_comp(9);
M_eqmin(10) = KOH*M_liq_comp(10);
M_eqmin(11) = KCO*M_liq_comp(11);

eqmin = M_eqmin.*MW;                % this is the bulk comosition of equilibrium magnesiowustite
eqmin = eqmin./(0.01*sum(eqmin));   % renormalized and in wt%

if eqmin(10) > 0.0075                  % prevent OH addition over saturation
    eqmin(10) = 0.0075;
    eqmin = eqmin./(0.01*sum(eqmin));
end

if eqmin(11) > 0.0001                  % prevent CO addition over saturation
    eqmin(11) = 0.0001;
    eqmin = eqmin./(0.01*sum(eqmin));
end
 
% now divide into MgO and FeO
 
Mgnum = M_eqmin(4)/(M_eqmin(4) + M_eqmin(3));
