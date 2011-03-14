function [eqmin, Mgnum] = spinel(liq_comp, MW)
% SPINEL calculates spinel compositional equilibrium
% (Mg, Fe)(Al)2(O)4
% liq_comp is in wt%
% MW is the molar weights of oxides
 
FeMgKd = 0.38;  % Kd = (Fe/Mg)min / (Fe/Mg)liq
 
M_liq_comp = liq_comp./MW;

KSm =0.001; 
KNd =0.001;
KTh =0.002;     
KU =0.012;     
KOH = 0.02;    % KOH = (OH)min/(OH)liq
KCO = 0.0005;

M_liq_comp = liq_comp./MW;
 
% use molar weights and Kds to calculate ratios of elements in mineral

FeMgliq = M_liq_comp(3)/M_liq_comp(4);
FeMgmin = FeMgKd*FeMgliq;

% Assum 1 Al2O3;

FirstFe = FeMgmin;
FirstMg = 1;
M_eqmin(1) = 0;
M_eqmin(2) = 1;
M_eqmin(3) = FirstFe/(FirstFe + FirstMg);
M_eqmin(4) = FirstMg/(FirstFe + FirstMg);
M_eqmin(5) = 0;
M_eqmin(6) = KSm*M_liq_comp(6);
M_eqmin(7) = KNd*M_liq_comp(7);
M_eqmin(8) = KTh*M_liq_comp(8);
M_eqmin(9) =  KU*M_liq_comp(9);
M_eqmin(10) = KOH*M_liq_comp(10);
M_eqmin(11) = KCO*M_liq_comp(11);

%M_eqmin
eqmin = M_eqmin.*MW;                % this is the bulk comosition of equilibrium garnet
eqmin = eqmin./(0.01*sum(eqmin));   % renormalized and in wt%

if eqmin(10) > 0.2                  % prevent OH addition over saturation
    eqmin(10) = 0.2;
    eqmin = eqmin./(0.01*sum(eqmin));
end

if eqmin(11) > 0.00001                  % prevent CO addition over saturation
    eqmin(11) = 0.00001;
    eqmin = eqmin./(0.01*sum(eqmin));
end

% now divide into (Ca)3(Al)2 garnet, (Mg,Fe)3(Al)2 garnet, [[(Mg, Fe)3 no Al garnet]]
% all with a basis of 3 silicas per unit, req. 2 Als per unit
 
Mgnum = M_eqmin(4)/(M_eqmin(4) + M_eqmin(3));
 