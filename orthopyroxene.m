function [eqmin, Mgnum, CaMg] = orthopyroxene(liq_comp, MW)

% ORTHOPYROXENE calculates equilibrium for orthopyroxene
% (Mg, Fe)2(Si)2(O)6
% liq_comp is in wt%
% MW is the molar weights of oxides
 
FeMgKd = 0.28;  % Kd = (Fe/Mg)min / (Fe/Mg)liq
CaMg = .02;       % Molar distrubution of Ca/Mg+Fe in pyroxene
 
M_liq_comp = liq_comp./MW;

KSm = 0.02;
KNd = 0.009; 
KTh = 0.13; 
KU  = 0.0004;
KOH = 0.02;    % KOH = (OH)min/(OH)liq 
KCO = 0.0004;

% use molar weights and Kds to calculate ratios of elements in mineral
 
FeMgliq = M_liq_comp(3)/M_liq_comp(4);
FeMgmin = FeMgKd*FeMgliq;
 
% Assume Fe = FeMgKd and Mg = 1 then normalize to 2(Fe,Mg) to 2 Si
 
FirstFe = FeMgmin;
FirstMg = 1;
 
M_eqmin(1) = 2;
M_eqmin(2) = 0;
M_eqmin(3) = 2*FirstFe/(FirstFe + FirstMg);
M_eqmin(4) = 2*FirstMg/(FirstFe + FirstMg);
M_eqmin(5) = (M_eqmin(4)+M_eqmin(3))*CaMg;
M_eqmin(6) = KSm*M_liq_comp(6);
M_eqmin(7) = KNd*M_liq_comp(7);
M_eqmin(8) = KTh*M_liq_comp(8);
M_eqmin(9) =  KU*M_liq_comp(9);
M_eqmin(10) = KOH*M_liq_comp(10);
M_eqmin(11) = KCO*M_liq_comp(11);
 
%M_eqmin
eqmin = M_eqmin.*MW;                % this is the bulk composition of equilibrium pyroxene
eqmin = eqmin./(0.01*sum(eqmin));   % renormalized and in wt%
%disp('pyroxene'); eqmin

if eqmin(10) > 0.1500                  % prevent OH addition over saturation
    eqmin(10) = 0.1500;
    eqmin = eqmin./(0.01*sum(eqmin));
end

if eqmin(11) > 0.00001                  % prevent CO addition over saturation
    eqmin(11) = 0.00001;
    eqmin = eqmin./(0.01*sum(eqmin));
end

Mgnum = M_eqmin(4)/(M_eqmin(4) + M_eqmin(3));

