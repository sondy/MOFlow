% Liquidprofile.m
% Takes liquid matrix from a profile calculation,
% along with P, and calcs density of liquid at each step
% Birch-Murnaghan solution for melt density with depth
% technique from Delano (1990)
% P is pressure in GPa
 
 
adjust = 2;    % amount to change rho if error is
   % too large, in kg/m^3
error = 0.05;   % max allowable P error in GPa, ends iterative loop
rhodP = 80; % d(rho)/dP in kg/m^3 per GPa, used for first guess
Kot = 22;    % bulk modulus for silicate liquids in GPa
Kpt = 6;     % pressure derivative of bulk modulus, varies from 5 to 7
 
 
for j = 1:i;
rhon = 1000*surfacemeltdensity(liquid(j,:), T(j));   % density of liquid at one atmosphere
rhonref = 1000*surfacemeltdensity(liquid(j,:), 1);
     ActPressure = P(j);
     rho = rhon + (ActPressure)*rhodP; % initial guess
CalcPressure = (3/2)*Kot*[(rho/rhon)^(7/3) - (rho/rhon)^(5/3)]*[1-...
        + (3/4)*(4 - Kpt)*[(rho/rhon)^(2/3) - 1]];
 
while abs(CalcPressure-ActPressure) >= error
     if (CalcPressure - ActPressure) > 0 % CalcPressure too high, rho too low
         rho = rho - adjust;
       elseif (CalcPressure - ActPressure) < 0
         rho = rho + adjust;
       end
CalcPressure = (3/2)*Kot*[(rho/rhon)^(7/3) - (rho/rhon)^(5/3)]*[1 - ...
        + (3/4)*(4 - Kpt)*[(rho/rhon)^(2/3) - 1]];
end
 
Dliquid(j) = rho;
Doliquid(j) = rhonref;
Dosolliquid(j) = rhon;
end
 
% Dliquid(j+1) = Dliquid(j) + (Dliquid(j) - Dliquid(j-1));
% Doliquid(j+1) = Doliquid(j) + (Doliquid(j) - Doliquid(j-1));
% Dosolliquid(j+1) = Dosolliquid(j) + (Dosolliquid(j) - Dosolliquid(j-1));
 
% figure(1);
% plot(Dliquid,r,'g');
%  
% figure(2);
% plot(Doliquid,r,'g');
%  
figure(3);
plot(Dosolliquid,r,'g');
 
figure(8);
plot(Dliquid, r, 'r');

% figure(25);
% plot(Dosolliquid, liqfrac)
% xlabel('liquid density at one atmosphere and liquidus temperature');
% ylabel('fraction that is still liquid');
% axis ij
% title(['density of evolving liquid at one atm for model:  ', name]);
 
figure(26);
plot(Doliquid, liqfrac)
xlabel('liquid density at one atmosphere and 1 degree');
ylabel('fraction that is still liquid');
axis ij
title(['density of evolving liquid at one atm and one degree for model:  ', name]);

% figure(27);
% plot(liquid,r);
% legend('SiO2','Al2O3','FeO','MgO','CaO','Sm','Nd','Lu','Hf');
% title('Mars moving solidus liquid composition');
% axis([0 60 1000 3500]);
 