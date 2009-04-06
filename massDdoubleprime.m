% massDdoubleprime.m
% 3/13/2009; Alessondra Springmann
% A routine to calculate the mass of the D'' layer via numerical
% integration

% Imagine a spherical shell with surface area SA = 4*pi*r^2
% A shell would have volume dV = 4*pi*r^2*dr and mass dM = rho(r)*dV
% The mass of this layer is therefore
% M = int_r1^r2 4*pi*rho(r)*r^2

% trapz is a numerical integration routine in Matlab

rho2_index = find(3317 < Dsolinv & Dsolinv < 3319); % density at the top of 
% the D'' layer

rho1_index = 1; % density at the CMB

radius1 = rinv(rho1_index); % radius at the CMB

radius2 = rinv(rho2_index); % radius at the top of the D'' layer

rhoDprimeprime = Dsolinv(rho1_index:rho2_index);

radiusDprimeprime = rinv(rho1_index:rho2_index);

massDprimeprime = trapz(radiusDprimeprime,...
    4*pi.*rhoDprimeprime.*radiusDprimeprime.^2);

fprintf('The mass of the D double prime layer is %2.3g kg. \n',...
    massDprimeprime)

Mearth = 5.9742e24;

percEarthMass = massDprimeprime*100/Mearth;

fprintf('This is %2.3g%% the mass of the Earth. \n', percEarthMass)

percMantleMass = massDprimeprime*100/Mantlemass;

fprintf('This is %2.3g%% the mass of the mantle. \n', percMantleMass)