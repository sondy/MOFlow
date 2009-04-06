function [eqmin, Mgnum] = olivine(liq_comp)

% OLIVINE calculates olivine equilibrium
% (Mg, Fe)2(Si)2(O)6
% liq_comp is in wt%
% MW is the molar weights of oxides
 
FeMgKd = 0.30;  % Kd = (Fe/Mg)min / (Fe/Mg)liq
 
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

KSm = 0.0006;               % from Snyder et al. 1992
KNd = 0.0001;
KU =  0.01;
KTh = 0.03;
KOH = 0.002;    % KOH = (OH)min/(OH)liq
KCO = 0.0005;

%KLu = 0.02;
%KHf = 0.01;

% use molar weights and Kds to calculate ratios of elements in mineral
 
FeMgliq = M_liq_comp(3)/M_liq_comp(4);
 
FeMgmin = FeMgKd*FeMgliq;

% Assume Fe = FeMgKd and Mg = 1 then normalize to 2(Fe,Mg) to 2 Si
 
FirstFe = FeMgmin;
FirstMg = 1;
M_eqmin(1) = 1;
M_eqmin(2) = 0;
M_eqmin(3) = 2*FirstFe/(FirstFe + FirstMg);
M_eqmin(4) = 2*FirstMg/(FirstFe + FirstMg);
M_eqmin(5) = 0;
M_eqmin(6) = KSm*M_liq_comp(6);
M_eqmin(7) = KNd*M_liq_comp(7);
M_eqmin(8) = KTh*M_liq_comp(8);
M_eqmin(9) =  KU*M_liq_comp(9);
M_eqmin(10) = KOH*M_liq_comp(10);
M_eqmin(11) = KCO*M_liq_comp(11);

% M_eqmin
eqmin = M_eqmin.*MW;                % this is the bulk composition of equilibrium olivine
eqmin = eqmin./(0.01*sum(eqmin));   % renormalized and in wt%

if eqmin(10) > 0.12                  % prevent OH addition over saturation
    eqmin(10) = 0.12;
    eqmin = eqmin./(0.01*sum(eqmin));
end

if eqmin(11) > 0.0005                 % prevent CO addition over saturation
    eqmin(11) = 0.0005;
    eqmin = eqmin./(0.01*sum(eqmin));
end
 
Mgnum = M_eqmin(4)/(M_eqmin(4) + M_eqmin(3)); 
% display(Mgnum)