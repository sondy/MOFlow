function [eqmin, Mgnum, Perc_Ca, Perc_AlMgFe, Perc_MgFe] = majorite(liq_comp, MW)
% MAJORITE calculates majorite equilibrium
% (Mg, Fe, Ca)3(Al)2(Si)3(O)12
% Order of elements = SiO2, Al2O3, FeO, MgO, CaO, Sm, Nd, Lu, Hf
% liq_comp is in wt%
% MW is the molar weights of oxides
 
% wt% aluminum and silica are from experiments from Tronnes and Frost (2002)
% otherwise determining alumina percent is not possible
 
FeMgKd = 0.45;  % Kd = (Fe/Mg)min / (Fe/Mg)liq
CaAlKd = 0.4;	% Kd = (Ca/Al)min / (Ca/Al)liq
SiMgKd = 1.2;	% Kd = (Si/Mg)min / (Si/Mg)liq
 
Alwt = 11; % 15 is experimental wt% from Bertka and Fei; Tronnes is different (4 and 51)
Siwt = 46; % experimental wt%
 
M_liq_comp = liq_comp./MW;

KSm = 0.3;
KNd = 0.05;
KTh = 0.02;
KU  = 0.02;
KOH = 0.003;    % KOH = (OH)min/(OH)liq from Bolfan-Casanova paper
KCO = 0.0005;

AlSi = (Alwt/MW(2))/(Siwt/MW(1));
M_liq_comp = liq_comp./MW;
 
% use molar weights and Kds to calculate ratios of elements in mineral
 
FeMgliq = M_liq_comp(3)/M_liq_comp(4);
CaAlliq = M_liq_comp(5)/(2*M_liq_comp(2));
SiMgliq = M_liq_comp(1)/M_liq_comp(4);
 
FeMgmin = FeMgKd*FeMgliq;
CaAlmin = CaAlKd*CaAlliq;
SiMgmin = SiMgKd*SiMgliq;
 
% assume 3 silicas per formula unit
% use Si/Mg to calculate moles of Mg
 
M_eqmin(1) = 3;
M_eqmin(4) = SiMgmin*M_eqmin(1);
 
% Fe follows from Mg
 
M_eqmin(3) = FeMgmin*M_eqmin(4);
 
% Calc Al from experimental ratios
 
 
M_eqmin(2) = AlSi*M_eqmin(1);
M_eqmin(5) = CaAlmin*(2*M_eqmin(2));
M_eqmin(6) = KSm*M_liq_comp(6);
M_eqmin(7) = KNd*M_liq_comp(7);
M_eqmin(8) = KTh*M_liq_comp(8);
M_eqmin(9) = KU*M_liq_comp(9);
M_eqmin(10) = KOH*M_liq_comp(10);
M_eqmin(11) = KCO*M_liq_comp(11);
 
eqmin = M_eqmin.*MW;                % this is the bulk comosition of 
                                    % equilibrium majorite
eqmin = eqmin./(0.01*sum(eqmin));   % renormalized and in wt%

if eqmin(10) > 0.0675                  % prevent OH addition over saturation
    eqmin(10) = 0.0675;
    eqmin = eqmin./(0.01*sum(eqmin));
end

if eqmin(11) > 0.0001                  % prevent CO addition over saturation
    eqmin(11) = 0.0001;
    eqmin = eqmin./(0.01*sum(eqmin));
end

% now divide into (Ca)3(Al)2 majorite, (Mg,Fe)3(Al)2 majorite, 
%     (Mg, Fe)3 no Al majorite
% all with a basis of 3 silicas per unit, req. 2 Als per unit
 
Mgnum = M_eqmin(4)/(M_eqmin(4) + M_eqmin(3));
if M_eqmin(5)/3 > M_eqmin(2)    % use up the aluminum in Ca majorite first
     TotCa = M_eqmin(2);         % if more Ca/3 than Al/2, use up all Al
     TotAlMgFe = 0;              % no Al left over for (Mg,Fe)3(Al)2 
     TotMgFe = (M_eqmin(3) + M_eqmin(4))/3;  % Fe and Mg go into (Mg, Fe)3 
else
     TotCa = M_eqmin(5)/3;       % otherwise use up all Ca
     TotAlMgFe = M_eqmin(2) - (M_eqmin(5)/3);    % and remainder of Al 
                                        % determines (Mg,Fe)3(Al)2 majorite
     TotMgFe = (M_eqmin(3) + M_eqmin(4))/3 - TotAlMgFe;   % the rest 
                    % of the Fe and Mg go into (M_eqmin(3) + M_eqmin(4))/3
end
Perc_Ca = TotCa/(TotCa + TotAlMgFe + TotMgFe);   % fraction of Ca, Al majorite
Perc_AlMgFe = TotAlMgFe/(TotCa + TotAlMgFe + TotMgFe);   % fraction 
                                            %     of Mg, Fe, Al majorite
Perc_MgFe = TotMgFe/(TotCa + TotAlMgFe + TotMgFe);   % fraction of 
                                            %  non-Al majorite
 
