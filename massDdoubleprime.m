% massDdoubleprime.m
% 3/13/2009; Alessondra Springmann
% A routine to calculate the mass of the D'' layer via numerical
% integration

% Imagine a spherical shell with surface area SA = 4*pi*r^2
% A shell would have volume dV = 4*pi*r^2*dr and mass dM = rho(r)*dV
% The mass of this layer is therefore
% M = int_r1^r2 4*pi*rho(r)*r^2

% trapz is a numerical integration routine in Matlab

rho2_index = find(3100 < Dsolinv & Dsolinv < 3200); % density at the top of 
% the D'' layer

rho1_index = 1; % density at the CMB

radius1 = rinv(rho1_index); % radius at the CMB

radius2 = rinv(50); % radius at the top of the D'' layer

% display(radius1)
% display(radius2)

rhoDprimeprime = mean(Dsolinv(rho1_index:rho2_index));

fprintf('The average density of the D" layer is %2.3g kg/m^3. \n',...
    rhoDprimeprime)

volumeDprimeprime = (4*pi/3)*(rinv(50).^3 - rinv(rho1_index).^3);

fprintf('The volume of the D" is %2.3g m^3. \n', volumeDprimeprime)

radiusDprimeprime = rinv(rho1_index:rho2_index);

massDprimeprime = trapz(radiusDprimeprime,...
    4*pi.*rhoDprimeprime.*radiusDprimeprime.^2);

fprintf('The mass of the D double prime layer is %2.3g kg. \n',...
    massDprimeprime)

percEarthMass = massDprimeprime*100/Mearth;

fprintf('This is %2.3g%% the mass of the Earth. \n', percEarthMass)

percMantleMass = massDprimeprime*100/mass_of_mantle;

fprintf('This is %2.3g%% the mass of the mantle. \n', percMantleMass)

%% Residual Liquid Density
densityResidual = new_mass_liquid./totalliquidvol;

%fprintf('The new mass liquid is %2.3g kg. \n', new_mass_liquid)

fprintf('The volume of the residual liquids is %2.3g kg. \n', totalliquidvol)

fprintf('The mass of the residual liquids is %2.3g kg. This is %2.3g%% the mass of the mantle. \n',...
    new_mass_liquid, (new_mass_liquid*100/mass_of_mantle))

fprintf('The density of the residual liquids is %2.3g kg/m^3 \n',...
    densityResidual)

%% Mass Balance - a test, of sorts

%fprintf('The mass of the mantle is %2.3g kg. \n', mass_of_mantle)

fprintf('The mass of the solidified ocean + residual liquids is %2.3g kg. \n',...
    (mass_solidified + new_mass_liquid))

silicateMass = (((R-DM).^3 - CMB.^3)/(R.^3 - CMB.^3))*mass_of_mantle;

% SO - solidified ocean; RL - residual liquids; US - unmelted silicates

fprintf('The mass of the SO, RL, and US is %3.4g%% of the present mantle mass. \n', ...
    ((silicateMass + mass_solidified + new_mass_liquid)*100)/(mass_of_mantle))

%% Soo, does liq_comp make sense?
% 12/4/2010

% 1: SiO2, 2: Al2O3, 3: FeO, 4: MgO, 5: CaO, 6: Sm, 7: Nd, 8: Th, 9: U, 10:
% OH, 11: C

%disp(liq_comp)

% Take the values from liq_comp and multiply them by the density of each
% species.  Does the resulting density make sense?

liq_comp_density = liq_comp.*mineral_density;

density_calc_res_liq = sum((liq_comp_density./100));

fprintf('The calculated density of the final layer is %2.3g kg/m^3. \n',...
    density_calc_res_liq)

mass_calc_res_liq = density_calc_res_liq*totalliquidvol;

% Well if it's close to 3600 kg/m^3, then yes.

%% Mass of the residual liquids combined with the \Dpp

mass_calc_residual = density_calc_res_liq * totalliquidvol;

mass_res_dpp = mass_calc_residual + massDprimeprime;