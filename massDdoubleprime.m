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