% massDdoubleprime.m
% 3/13/2009; Alessondra Springmann
% A routine to calculate the mass of the D'' layer via numerical
% integration

% Imagine a spherical shell with surface area SA = 4*pi*r^2
% A shell would have volume dV = 4*pi*r^2*dr and mass dM = rho(r)*dV
% The mass of this layer is therefore
% M = int_r1^r2 4*pi*rho(r)*r^2

% trapz is a numerical integration routine in Matlab

% February 2011: LTET altered to add last liquid layers into both mass and
% volume of D" candidate; assumes all last liquids sink, and assumes they
% have the density of the average of the solid D"

rho2_index = find(3317 < Dsolinv & Dsolinv < 3320); % density at the top of D" - chosen to find kink in density
rho2_index = rho2_index(1);
rho1_index = 1; % density at the CMB

radius1 = rinv(rho1_index); % radius at the CMB
radius2 = rinv(rho2_index); % radius at the top of the D'' layer

rhoDprimeprime = Dsolinv(rho1_index:rho2_index);    % vector of D" densities
avgrhoD = mean(rhoDprimeprime);                  % average D" density
radiusDprimeprime = rinv(rho1_index:rho2_index);    % vector of D" radii

totalDvol = (4/3)*pi*(radius2^3 - radius1^3); % total D" volume
totalliquidvol = (4/3)*pi*(R^3 - r(maxstep)^3); % total liquid unsolidified at top of MO

massDprimeprime = trapz(radiusDprimeprime,...
    4*pi.*rhoDprimeprime.*radiusDprimeprime.^2) + ...
    totalliquidvol*avgrhoD;

fprintf('The mass of the D double prime layer is %2.3g kg. \n',...
    massDprimeprime)

Mearth = 5.9742e24;
Mantlemass = (R^3)/(R^3 - CMB^3)*4.032e+024;

percEarthMass = massDprimeprime*100/Mearth;

fprintf('This is %2.3g%% the mass of the Earth. \n', percEarthMass)

percMantleMass = massDprimeprime*100/Mantlemass;

fprintf('This is %2.3g%% the mass of the mantle. \n', percMantleMass)

% calculates average Sm, Nd, U, Th content for EER (=D") layers of "max" steps above CMB
max = rho2_index;

for k = 2:max
    delavgSm(k) = solidinv(k,6)*(4/3)*pi*(rinv(k)^3 - rinv(k-1)^3)/totalDvol;
    delavgNd(k) = solidinv(k,7)*(4/3)*pi*(rinv(k)^3 - rinv(k-1)^3)/totalDvol;
    delavgU(k) = solidinv(k,9)*(4/3)*pi*(rinv(k)^3 - rinv(k-1)^3)/totalDvol;
    delavgTh(k) = solidinv(k,8)*(4/3)*pi*(rinv(k)^3 - rinv(k-1)^3)/totalDvol;
end

avgNdEER = sum(delavgNd);   % for insertion into Rick Carlson's spreadsheet
avgSmEER = sum(delavgSm);

% disp(['D" Nd wt% is ', num2str(avgNdEER),' and Sm wt% is ', num2str(avgSmEER)])
% % these two include all the residual liquid from the surface
avgNdEERwLiq = (totalliquidvol/(totalliquidvol+totalDvol))*liquid(maxstep,7) + (totalDvol/(totalliquidvol+totalDvol))*avgNdEER;
avgSmEERwLiq = (totalliquidvol/(totalliquidvol+totalDvol))*liquid(maxstep,6) + (totalDvol/(totalliquidvol+totalDvol))*avgSmEER;
% disp(['With all final liquids D" Nd wt% is ', num2str(avgNdEERwLiq),' and Sm wt% is ', num2str(avgSmEERwLiq)])
% 
avgU = sum(delavgU);
avgTh = sum(delavgTh);
EERUfracoftotal = avgU*totalDvol/(liquid(1,9)*Mantlevolume);  % for comparison with Carlson's estimates
EERThfracoftotal = avgTh*totalDvol/(liquid(1,8)*Mantlevolume);    % of what U and Th fraction must be in D"
% disp(['D" fraction of total Earth U is ', num2str(EERUfracoftotal),' and Th fraction is ', num2str(EERThfracoftotal)])
% 
% %% with liquids
EERUfracoftotalwliq = (avgU*totalDvol + liquid(maxstep,9)*totalliquidvol)/(liquid(1,9)*Mantlevolume);
EERThfracoftotalwliq = (avgTh*totalDvol + liquid(maxstep,8)*totalliquidvol)/(liquid(1,8)*Mantlevolume);
% disp(['With all final liquids D" fraction of total Earth U is ', num2str(EERUfracoftotalwliq),' and Th fraction is ', num2str(EERThfracoftotalwliq)])
% 
avgNdEDR = mean(solidinv(max:990,7));
avgSmEDR = mean(solidinv(max:990,6));
% disp(['Mantle minus D" Nd wt% is ', num2str(avgNdEDR),' and Sm wt% is ', num2str(avgSmEDR)])