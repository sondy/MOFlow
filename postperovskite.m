function [eqmin, Mgnum] = postperovskite(liq_comp)

% POSTPEROVSKITE calculates equilibrium compositions
% (Mg, Fe)(Si)(O)3
% liq_comp is in wt%
% MW is the molar weights of oxides
 
FeMgKd = 0.4;  % Kd = (Fe/Mg)min / (Fe/Mg)liq
 
M_liq_comp = liq_comp./MW;

load kd_postperovskite.dat
kd_pp = kd_postperovskite;

KSm = kd_pp(1);%0.15;
KNd = kd_pp(2);%0.05; 
KTh = kd_pp(3);%0.5; 
KU  = kd_pp(4);%0.5;
KOH = kd_pp(5);%0.0001;    % KOH = (OH)min/(OH)liq 
KCO = kd_pp(6);%0.0005;

% use molar weights and Kds to calculate ratios of elements in mineral
 
FeMgliq = M_liq_comp(3)/M_liq_comp(4);
FeMgmin = FeMgKd*FeMgliq;
 
% Assume Fe = FeMgKd and Mg = 1 then normalize to 1 (Fe,Mg) to 1 Si
 
FirstFe = FeMgmin;
FirstMg = 1;
 
M_eqmin(1) = 1;
M_eqmin(2) = 0;
M_eqmin(3) = FirstFe/(FirstFe + FirstMg);
M_eqmin(4) = FirstMg/(FirstFe + FirstMg);
M_eqmin(5) = 0;
M_eqmin(6) = KSm*M_liq_comp(6);
M_eqmin(7) = KNd*M_liq_comp(7);
M_eqmin(8) = KTh*M_liq_comp(8);
M_eqmin(9) = KU*M_liq_comp(9);
M_eqmin(10) = KOH*M_liq_comp(10);
M_eqmin(11) = KCO*M_liq_comp(11);
 
%M_eqmin
eqmin = M_eqmin.*MW;                % this is the bulk composition of equilibrium pyroxene
eqmin = eqmin./(0.01*sum(eqmin));   % renormalized and in wt%

if eqmin(10) > 0.01                  % prevent OH addition over saturation
    eqmin(10) = 0.01;
    eqmin = eqmin./(0.01*sum(eqmin));
end

if eqmin(11) > 0.0001                  % prevent CO addition over saturation
    eqmin(11) = 0.0001;
    eqmin = eqmin./(0.01*sum(eqmin));
end

Mgnum = M_eqmin(4)/(M_eqmin(4) + M_eqmin(3));

