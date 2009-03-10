function [eqmin, Mgnum, Perc_AlCa, Perc_AlMgFe] = garnet(liq_comp)
% GARNET calculates garnet equilibrium
% (Mg, Fe, Ca)3(Al)2(Si)3(O)12
% liq_comp is in wt%
% MW is the molar weights of oxides
 
% assume 2 Als per formula unit no matter what, also 4 wt% CaO no matter  what
 
FeMgKd = 0.48;  % Kd = (Fe/Mg)min / (Fe/Mg)liq
CaAlKd = 0.4;
SiMgKd = 0.8;
 
% Alwt = 22; % experimental wt% from Bertka and Fei
% Siwt = 41; % experimental wt%
 
MW = [  60.09,...   % SiO2   (1)   g/mol for oxides
       101.96,...   % Al2O3  (2)
        71.846,...  % FeO    (3)
        40.311,...  % MgO    (4)
        56.077,...  % CaO    (5)
       150.36,...   % Sm     (6) (these in molar ppm after division)
       144.24,...   % Nd     (7)
       232.0381,... % Th     (8)
       238.0289,... % U      (9)
       17.00,...    % OH     (10)   g/mol again
       28.00];      % CO     (11)
 
       %174.967,...  % Lu     (8)
       %178.49,...   % Hf     (9)   
 
M_liq_comp = liq_comp./MW;

KSm =3;            
KNd =0.05;
KTh =0.003;   
KU =0.01;    
KOH = 0.0008;    % KOH = (OH)min/(OH)liq
KCO = 0.0002;

% AlSi = (Alwt/MW(2))/(Siwt/MW(1));
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
Mg = SiMgmin*M_eqmin(1);
 
% Fe follows from Mg
 
Fe = FeMgmin*Mg;
 
% Assum 2 Als (1 Al2O3); assume 4% CaO (experimental)
 
M_eqmin(2) = 1;
M_eqmin(5) = 0.30;
% M_eqmin(5) = CaAlmin*(M_eqmin(2));
% normalize Fe, Mg
M_eqmin(3) = 3*Fe/(Fe+Mg);
M_eqmin(4) = 3*Mg/(Fe+Mg);
M_eqmin(6) = KSm*M_liq_comp(6);
M_eqmin(7) = KNd*M_liq_comp(7);
M_eqmin(8) = KTh*M_liq_comp(8);
M_eqmin(9) =  KU*M_liq_comp(9);
M_eqmin(10) = KOH*M_liq_comp(10);
M_eqmin(11) = KCO*M_liq_comp(11);

%M_eqmin
eqmin = M_eqmin.*MW;                % this is the bulk comosition of equilibrium garnet
eqmin = eqmin./(0.01*sum(eqmin));   % renormalized and in wt%

if eqmin(10) > 0.0700                  % prevent OH addition over saturation
    eqmin(10) = 0.0700;
    eqmin = eqmin./(0.01*sum(eqmin));
end

if eqmin(11) > 0.00015                  % prevent CO addition over saturation
    eqmin(11) = 0.00015;
    eqmin = eqmin./(0.01*sum(eqmin));
end

% now divide into (Ca)3(Al)2 garnet, (Mg,Fe)3(Al)2 garnet, [[(Mg, Fe)3 no Al garnet]]
% all with a basis of 3 silicas per unit, req. 2 Als per unit
 
Mgnum = M_eqmin(4)/(M_eqmin(4) + M_eqmin(3));
TotAlCa = M_eqmin(5);
TotAlMgFe = (M_eqmin(3) + M_eqmin(4));
 
% if M_eqmin(5)/3 > M_eqmin(2)    % use up the aluminum in Ca garnet first
%      TotCa = M_eqmin(2);         % if there is more Ca/3 than Al/2, use up all Al
%      TotAlMgFe = 0;              % there is no Al left over for (Mg,Fe)3(Al)2 garnet
%      TotMgFe = (M_eqmin(3) + M_eqmin(4))/3;  % all Fe and Mg go into (Mg, Fe)3 no Al garnet
% else
%      TotCa = M_eqmin(5)/3;       % otherwise use up all Ca
%      TotAlMgFe = M_eqmin(2) - (M_eqmin(5)/3);    % and remainder of Al determines (Mg,Fe)3(Al)2 garnet
%      TotMgFe = (M_eqmin(3) + M_eqmin(4))/3 - TotAlMgFe;   % the rest of the Fe and Mg go into (M_eqmin(3) + M_eqmin(4))/3
% end
Perc_AlCa = TotAlCa/(TotAlCa + TotAlMgFe);   % fraction of Ca, Al garnet
Perc_AlMgFe = TotAlMgFe/(TotAlCa + TotAlMgFe);   % fraction of Mg, Fe, Al garnet
%Perc_MgFe = TotMgFe/(TotCa + TotAlMgFe + TotMgFe);   % fraction of non-Al garnet
 