function [eqmin, Mgnum] = gammamin(liq_comp)
% GAMMA (GAMMA OLIVINE, RINGWOODITE) calculates gamma equilibrium
% (Mg, Fe)2(Si)1(O)4
% given liquid composition with five components
% Order of elements = SiO2, Al2O3, FeO, MgO, CaO
% liq_comp is in wt%
% MW is the molar weights of oxides
 
FeMgKd = 0.7;  % Kd = (Fe/Mg)min / (Fe/Mg)liq
 
M_liq_comp = liq_comp./MW;

KSm = 0.0006;
KNd = 0.0001;
KTh = 0.02; % should be changed
KU  = 0.01; % who knows?
KOH = 0.03;    % KOH = (OH)min/(OH)liq from Bolfan-Casanova paper
KCO = 0.0005;

% use molar weights and Kds to calculate ratios of elements in mineral
 
FeMgliq = M_liq_comp(3)/M_liq_comp(4);
FeMgmin = FeMgKd*FeMgliq;
 
% Assume Fe = FeMgKd and Mg = 1 then normalize to 2(Fe,Mg) to 1 Si
 
FirstFe = FeMgmin;
FirstMg = 1;
M_eqmin(1) = 1;
M_eqmin(3) = 2*FirstFe/(FirstFe + FirstMg);
M_eqmin(4) = 2*FirstMg/(FirstFe + FirstMg);
M_eqmin(5) = 0;
M_eqmin(2) = 0;
M_eqmin(6) = KSm*M_liq_comp(6);
M_eqmin(7) = KNd*M_liq_comp(7);
M_eqmin(8) = KTh*M_liq_comp(8);
M_eqmin(9) = KU*M_liq_comp(9);
M_eqmin(10) = KOH*M_liq_comp(10);
M_eqmin(11) = KCO*M_liq_comp(11);
 
%M_eqmin
eqmin = M_eqmin.*MW;                % this is the bulk comosition of equilibrium spinel
eqmin = eqmin./(0.01*sum(eqmin));   % renormalized and in wt%

if eqmin(10) > 2.5                  % prevent OH addition over saturation
    eqmin(10) = 2.5;
    eqmin = eqmin./(0.01*sum(eqmin));
end

if eqmin(11) > 0.000003                  % prevent CO addition over saturation
    eqmin(11) = 0.000003;
    eqmin = eqmin./(0.01*sum(eqmin));
end
 
Mgnum = M_eqmin(4)/(M_eqmin(4) + M_eqmin(3)); 
