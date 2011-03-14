function [eqmin, An] = plagioclase(liq_comp)
% PLAGIOCLASE calculates equilibrium for plagioclase
% Ca(Al)2(Si)2(O)6 (anorthite) and NaAlSi3O8 (albite)
% given liquid composition with seven components
% Order of elements = SiO2, Al2O3, FeO, MgO, CaO, TiO2, Cr2O3
% liq_comp is in wt%
% MW is the molar weights of oxides
 
An = 1;       % Fraction of plag that is anorthite; set to one as long as
                % there is no sodium in the calculation
 
M_liq_comp = liq_comp./MW;

KSm = 0.03;              
KNd = 0.04;
KTh = 0.12;
KU  = 0.01;
KOH = 0.001;    % KOH = (OH)min/(OH)liq 
KCO = 0.0005;

% Make up plag by adding anorthite and albite
  
M_eqmin(1) = An*2 + (1-An)*3;
M_eqmin(2) = An*1 + (1-An)*0.5;
M_eqmin(3) = 0;
M_eqmin(4) = 0;
M_eqmin(5) = An*1 + (1-An)*0;
M_eqmin(6) = KSm*M_liq_comp(6);
M_eqmin(7) = KNd*M_liq_comp(7);
M_eqmin(8) = KTh*M_liq_comp(8);
M_eqmin(9) =  KU*M_liq_comp(9);
M_eqmin(10) = KOH*M_liq_comp(10);
M_eqmin(11) = KCO*M_liq_comp(11);
 
%M_eqmin
eqmin = M_eqmin.*MW;                % this is the bulk composition of equilibrium plag
eqmin = eqmin./(0.01*sum(eqmin));   % renormalized and in wt%
%disp('plagioclase'); eqmin

if eqmin(10) > 0.0510                 % prevent OH addition over saturation
    eqmin(10) = 0.0510;
    eqmin = eqmin./(0.01*sum(eqmin));
end

if eqmin(11) > 0.0001                  % prevent CO addition over saturation
    eqmin(11) = 0.0001;
    eqmin = eqmin./(0.01*sum(eqmin));
end
