function [eqmin, Mgnum, Perc_Al, Perc_Ca, Perc_MgFe] = perovskite(liq_comp)
% PEROVSKITE calculates equilibrium perovskite
% (Mg, Fe, Ca)3(Al)2(Si)3(O)12
% given liquid composition with five components
% liq_comp is in wt%
% MW is the molar weights of oxides
 
FeMgKd = 0.46;  % Kd = (Fe/Mg)min / (Fe/Mg)liq
CaAlKd = 0.5;
SiMgKd = 1.1;
 
Alwt = 5; %mol% of Al2O3; experimental Al2O3 wt% from Bertka and Fei: 4; Tronnes is similar; Panero (2006) higher: 15 mol%
Siwt = 52; % experimental SiO2 wt%
 
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

AlSi = (Alwt/MW(2))/(Siwt/MW(1));
M_liq_comp = liq_comp./MW;

KSm = 0.15; 
KNd = 0.05;
KTh = 7; % corgne & wood
KU =  7;
KOH = 0.0001;    % KOH = (OH)min/(OH)liq
KCO = 0.0005; 

% use molar weights and Kds to calculate ratios of elements in mineral
 
FeMgliq = M_liq_comp(3)/M_liq_comp(4);
CaAlliq = M_liq_comp(5)/(M_liq_comp(2));
SiMgliq = M_liq_comp(1)/M_liq_comp(4);
 
FeMgmin = FeMgKd*FeMgliq;
CaAlmin = CaAlKd*CaAlliq;
SiMgmin = SiMgKd*SiMgliq;

% assume 3 silicas per formula unit
% use Si/Mg to calculate moles of Mg
 
M_eqmin(1) = 3;
M_eqmin(4) = M_eqmin(1)/SiMgmin;
 
% Fe follows from Mg
 
M_eqmin(3) = FeMgmin*M_eqmin(4);
 
% Calc Al from experimental ratios compared to Si
 
M_eqmin(2) = AlSi*M_eqmin(1);   % moles of Al2O3
M_eqmin(5) = CaAlmin*2*M_eqmin(2);
M_eqmin(6) = KSm*M_liq_comp(6);
M_eqmin(7) = KNd*M_liq_comp(7);
M_eqmin(8) = KTh*M_liq_comp(8);
M_eqmin(9) = KU*M_liq_comp(9);
M_eqmin(10) = KOH*M_liq_comp(10);
M_eqmin(11) = KCO*M_liq_comp(11);

eqmin = M_eqmin.*MW;                % this is the bulk composition of equilibrium perovskite
eqmin = eqmin./(0.01*sum(eqmin));  % renormalized and in wt%

OHsat = (eqmin(3)*0.0015 + eqmin(4)*0.0010 + eqmin(5)*0.004)/[eqmin(3)+eqmin(4)+eqmin(5)]; % saturation depends strongly on Fe, Mg, Ca
if eqmin(10) > OHsat               % prevent OH addition over saturation
    eqmin(10) = OHsat;                
    eqmin = eqmin./(0.01*sum(eqmin));
end

if eqmin(11) > 0.000003                  % prevent CO addition over saturation
    eqmin(11) = 0.000003;
    eqmin = eqmin./(0.01*sum(eqmin));
end

% now divide into (Mg, Fe)3(Al)2 perovskite, [CaSiO3]3 Ca perovskite, and [(Mg,Fe)SiO3]3 non-Al, non-Ca perovskite
% req. 2 Als or 3 Cas per unit
 
Mgnum = M_eqmin(4)/(M_eqmin(4) + M_eqmin(3));
TotMgFe = (3)*M_eqmin(2);        % moles of Mg+Fe that go into Al perovskite to charge balance Al
RemMgFe = (M_eqmin(3) + M_eqmin(4)) - TotMgFe;  % moles of Mg+Fe left to make non-Al, non-Ca perovskite
Perc_Al = (M_eqmin(2))/((M_eqmin(2)) + (RemMgFe) + (M_eqmin(5)));   % fraction of Al perovskite
Perc_Ca = (M_eqmin(5))/((M_eqmin(2)) + (RemMgFe) + (M_eqmin(5)));   % fraction of Ca perovskite
Perc_MgFe = (RemMgFe)/((M_eqmin(2)) + (RemMgFe) + (M_eqmin(5)));   % fraction of MgFe perovskite
 

