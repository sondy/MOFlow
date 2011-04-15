function [eqmin, Mgnum, Perc_Al, Perc_Ca, Perc_MgFe] = perovskite(liq_comp, MW, Kd_p_Ca, Kd_p_MgFe)
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

AlSi = (Alwt/MW(2))/(Siwt/MW(1));
M_liq_comp = liq_comp./MW;

load kd_perovskite.dat
kd_p = kd_perovskite;

% KSm = 0.15; 
% KNd = 0.05;

% KSm_Ca   = Kd_p_Ca(1); %kd_p(1,1);%9;
% KSm_MgFe = Kd_p_MgFe(1); %kd_p(1,2);%0.05;
% KSm_Al = KSm_MgFe;
% 
% KNd_Ca   = Kd_p_Ca(2); %kd_p(2,1);%7;
% KNd_MgFe = Kd_p_MgFe(2); %kd_p(2,2);%0.016;
% KNd_Al = KNd_MgFe;
% 
% KTh_Ca = Kd_p_Ca(3); %kd_p(3,1);%10;
% KTh_MgFe = Kd_p_MgFe(3); %kd_p(3,2);%0.005;
% KTh_Al = KTh_MgFe;
% 
% KU_Ca = Kd_p_Ca(4); %kd_p(4,1);%8;
% KU_MgFe = Kd_p_MgFe(4); %kd_p(4,2);%0.025;
% KU_Al = KU_MgFe;

KSm_Ca   = kd_p(1,1);%9;
KSm_MgFe = kd_p(1,2);%0.05;
KSm_Al = KSm_MgFe;

KNd_Ca   = kd_p(2,1);%7;
KNd_MgFe = kd_p(2,2);%0.016;
KNd_Al = KNd_MgFe;

KTh_Ca = kd_p(3,1);%10;
KTh_MgFe = kd_p(3,2);%0.005;
KTh_Al = KTh_MgFe;

KU_Ca = kd_p(4,1);%8;
KU_MgFe = kd_p(4,2);%0.025;
KU_Al = KU_MgFe;

% KTh = 0.01; % corgne et al. 2004
% KU =  0.03; % " "

% KOH_Ca = Kd_p_Ca(5);
% KOH_MgFe = Kd_p_MgFe(5);
% 
% KCO_Ca = Kd_p_Ca(6);
% KCO_MgFe = Kd_p_MgFe(6);

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

Mgnum = M_eqmin(4)/(M_eqmin(4) + M_eqmin(3));
TotMgFe = (3)*M_eqmin(2);        % moles of Mg+Fe that go into Al perovskite to charge balance Al
RemMgFe = (M_eqmin(3) + M_eqmin(4)) - TotMgFe;  % moles of Mg+Fe left to make non-Al, non-Ca perovskite
Perc_Al = (M_eqmin(2))/((M_eqmin(2)) + (RemMgFe) + (M_eqmin(5)));   % fraction of Al perovskite
Perc_Ca = (M_eqmin(5))/((M_eqmin(2)) + (RemMgFe) + (M_eqmin(5)));   % fraction of Ca perovskite
Perc_MgFe = (RemMgFe)/((M_eqmin(2)) + (RemMgFe) + (M_eqmin(5)));   % fraction of MgFe perovskite

% Perc_Al = 0.10;
% Perc_Ca = 0.02;
% Perc_MgFe = 0.88;

% M_eqmin(6) = KSm*M_liq_comp(6); % Sm
% M_eqmin(7) = KNd*M_liq_comp(7); % Nd

M_eqmin(6) = M_liq_comp(6)*(Perc_Al*KSm_Al + Perc_Ca*KSm_Ca +...
    Perc_MgFe*KSm_MgFe);
M_eqmin(7) = M_liq_comp(7)*(Perc_Al*KNd_Al + Perc_Ca*KNd_Ca +...
    Perc_MgFe*KNd_MgFe);

M_eqmin(8) = M_liq_comp(8)*(Perc_Al*KTh_Al + Perc_Ca*KTh_Ca +...
    Perc_MgFe*KTh_MgFe);
M_eqmin(9) = M_liq_comp(9)*(Perc_Al*KU_Al + Perc_Ca*KU_Ca +...
    Perc_MgFe*KU_MgFe);

% M_eqmin(8) = KTh*M_liq_comp(8);
% M_eqmin(9) = KU*M_liq_comp(9);

% KOH = Perc_MgFe*KOH_MgFe + Perc_Ca*KOH_Ca;
% KCO = Perc_MgFe*KCO_MgFe + Perc_Ca*KCO_Ca;

M_eqmin(10) = KOH*M_liq_comp(10);
M_eqmin(11) = KCO*M_liq_comp(11);

eqmin = M_eqmin.*MW;                % this is the bulk composition of equilibrium perovskite
eqmin = eqmin./(0.01*sum(eqmin));  % renormalized and in wt%

OHsat = (eqmin(3)*0.0015 + eqmin(4)*0.0010 + eqmin(5)*0.004)./(eqmin(3)+eqmin(4)+eqmin(5)); % saturation depends strongly on Fe, Mg, Ca
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