function [Dliquid, Doliquid, Dosolliquid]=NewBirchMurnliquid(intliq, T, P)
% NewBirchMurnliquid.m
% Takes liquid composition and T and P
% along with P, and calcs density of liquid
% Birch-Murnaghan solution for melt density with depth
% technique from Delano (1990)
% P is pressure in GPa
 
%display('In NewBirchMurnliquid')

% adjust = 1;    % amount to change rho if error is
   % too large, in kg/m^3
% error = 0.05;   % max allowable P error in GPa, ends iterative loop
rhodP = 100; % d(rho)/dP in kg/m^3 per GPa, used for first guess
Kot = 22;    % bulk modulus for silicate liquids in GPa
Kpt = 6;     % pressure derivative of bulk modulus, varies from 5 to 7

rhon = 1000*surfacemeltdensity(intliq, T);
rhonref = 1000*surfacemeltdensity(intliq, 1);

ActPressure = P;
rho = rhon + (ActPressure)*rhodP; % initial guess
     
% CalcPressure = (3/2)*Kot*((rho/rhon)^(7/3) - (rho/rhon)^(5/3))*(1-...
%         + (3/4)*(4 - Kpt)*((rho/rhon)^(2/3) - 1));

% rho_ratio = rho/rhon
% 

if imag(Kot) ~= 0 || imag(Kpt) ~= 0 || imag(ActPressure) ~= 0 || imag(P) ~= 0 
    display('Something is imaginary in NewBirchMurnliquid.  Oh dear.')
end

Pfunc = @(rho_ratio)(3/2).*Kot.*((rho_ratio).^(7/3) - ...
    (rho_ratio).^(5/3)).*(1 - (3/4).*(4 - Kpt).*((rho_ratio).^(2/3) - 1))...
    - ActPressure;     

% rho/rhon

rhon = rho/(fzero(Pfunc,rho/rhon));
 
% while abs(CalcPressure-ActPressure) >= error
%      if (CalcPressure - ActPressure) > 0 % CalcPressure too high, rho too low
%          rho = rho - adjust;
%      elseif (CalcPressure - ActPressure) < 0
%          rho = rho + adjust;
%      end
% CalcPressure = (3/2)*Kot*((rho/rhon)^(7/3) - (rho/rhon)^(5/3))*(1 - ...
%         + (3/4)*(4 - Kpt)*((rho/rhon)^(2/3) - 1));
% end 

Dliquid = rho;
Doliquid = rhonref;
Dosolliquid = rhon;


 
% Dliquid(j+1) = Dliquid(j) + (Dliquid(j) - Dliquid(j-1));
% Doliquid(j+1) = Doliquid(j) + (Doliquid(j) - Doliquid(j-1));
% Dosolliquid(j+1) = Dosolliquid(j) + (Dosolliquid(j) - Dosolliquid(j-1));
 
% figure(1);
% plot(Dliquid,r,'g');
%  
% figure(2);
% plot(Doliquid,r,'g');
%  
% figure(3);
% plot(Dosolliquid,r,'g');
 
% figure(8);
% plot(Dliquid, r, 'r');
%  
