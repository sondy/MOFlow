function rhon = surfacemeltdensity(liq_comp, Tpot)
% calculates density of a silicate melt as a function of
% composition, temperature, and pressure.  Composition
% and temperature dependence are from Lange and Carmichael (1987).
% Pressure dependence, which includes only the first derivative of
% V with respect to T, is from Kress and Carmichael (1991).
% CO2 values from Liu and Lange (2003) CMP 146 p370
% H2O from Ochs and Lange 1997, 1999
% The expression is believed to be accurate for pressures
% up to 2 - 3 GPa.
% Used in magma ocean calculations only for P = 1 atm; Birch-M. after.
 
wtpct(1) = liq_comp(1);     % SiO2
wtpct(2) = 0;               % TiO2
wtpct(3) = liq_comp(2);     % Al2O3
wtpct(4) = liq_comp(3);     % FeO
wtpct(5) = liq_comp(4);     % MgO
wtpct(6) = liq_comp(5);     % CaO
wtpct(7) = 0;               % Na2O
wtpct(8) = 0;               % K2O
wtpct(9) = 0;               % Na2O - TiO2
wtpct(10) = -liq_comp(2);   % Na2O - Al2O3
wtpct(11) = liq_comp(10);   % H2O
wtpct(12) = liq_comp(11);   % CO2
 
T = Tpot + 298;             % Kelvins

%    SiO2    TiO2     Al2O3,   FeO,    MgO,    CaO   Na2O    K2O  Na2O-TiO2 Na2O-Al2O3  H2O     CO2
MW = [60.09, 79.879, 101.96, 71.846, 40.304, 56.077, 61.979, 94.196, 1,          1,       18.99, 44.00];
mol = wtpct./MW;
X = mol/sum(mol);

%       SiO2    TiO2   Al2O3, FeO,    MgO,  CaO   Na2O    K2O  Na2O-TiO2 Na2O-Al2O3  H2O     CO2
V1673 = [26.90, 23.16, 37.11, 13.65, 11.45, 16.57, 28.78, 45.84, 20.21,    0,       27.7, 28.00];
% partial molar volume [cm3/mol]
 
dVdT = [0, 7.24e-3, 2.63e-3, 2.92e-3, 2.62e-3, 2.92e-3, 7.41e-3, 11.91e-3, 0, 0, 9.5e-3, 4e-3];
dVdP1673 = [-1.89e-4, -2.31e-4, -2.26e-4, -0.45e-4, 0.27e-4, 0.34e-4, -2.4e-4, -6.75e-4, 0, 10.18e-4, -3.2e-4, -3e-4];
dVdPdT = [1.3e-7, 0, 2.7e-7, -1.8e-7, -1.3e-7, -2.9e-7, -6.6e-7, -14.5e-7, 0, 0, 0, -4e-7];
 
mass = MW*X';   % mass of 1 mol melt
Tref = 1673;    % ref. temperature to which partial molar quantities are referred - don't change!
 
Pmin = 0.5; %0;       % bars
Pmax = 1.5; %1000;    % should be good to Pmax of 30,000
numsteps = 2; %10;
 
P = [Pmin:(Pmax-Pmin)/(numsteps-1):Pmax];
 
for i = 1:numsteps
     Vpart(i,:) = V1673 + dVdT*(T-Tref) + (dVdP1673+dVdPdT*(T-Tref))*P(i);   % partial molar vols of components
     V(i) = Vpart(i,:)*X';   % molar volume (cm3)
     rho(i) = mass/V(i);     % density 9g/cm3)
end
 
rhon = rho(1);  % surface liquid density to pass to BirchMurnliquid, in kg/m3
 