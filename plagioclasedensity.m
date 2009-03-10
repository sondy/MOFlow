function density = plagioclasedensity(An, P, T)
 
% calculates density of plagioclase
% calculate Vo from Vo298 and integral for alpha dT
% calculate V(P,T) from Birch-Murnaghan equation at given P (molar volumes at P and T)
% using KT and KpT, each calculated from their own derivatives
% calculate density from V(P,T) and molar weights of each variety of pyroxene
% physical parameters from Holland and Powell (1998) 
 
% at first, physical parameter matrix is size (7 x 3)
% where each row is a parameter
% in order 1 = Vo298, 2 = KT298, 3 = dK/dP, 4 = dK/dT, 5 = ao, 6 - a1, 7 = a2 
% below parameters are added: 8 = Vo(T), 9 = KT, 10 = V(P,T) in m3/mole
%                               11 = kg/mol of formula units of phase, 12 = density of phase in kg/m3
%                               13 = calculated alpha, /deg

% physical parameters for Ca(Al)2(Si)2(O)6 (anorthite) and NaAlSi3O8 (albite)

par(1,1) = An*(101) + (1-An)*(100); % Vo298, cm3/mole 
par(2,1) = An*(91.9) + (1-An)*59.3;      % KT298, GPa
par(3,1) = 4;           % dK/dP this is the default value used by Holland and Powell
par(4,1) = 1e-4;        % dK/dT this is the default value used by Holland and Powell
par(5,1) = An*(0.0000238) + (1-An)*(0.0000456);   % ao
% par(6,1) = 2.694e-9;               % a1
% par(7,1) = -0.5588;               % a2
par(11,1) = ((An*1 + (1-An)*0)*56.077 + (An*1 + (1-An)*1)*101.96 + (An*2 + (1-An)*3)*60.09)/1000;  % kg/moles of formula
            %   Ca                          Al                          Si
 
%    par(13,i) = (par(5,i)*(T+298) + (par(6,i)/2)*((T+298)^2) - par(7,i)*(1/(T+298)))...
%         - (par(5,i)*(298) + (par(6,i)/2)*(298^2) - par(7,i)*(1/298)); % thermal expansivity
%   par(13,1) = (par(5,1)); % thermal expansivity
   par(13,1) = par(5,1)*(T+298) - par(5,1)*(298); % thermal expansivity
   par(8,1) = par(1,1)*exp(par(13,1));        % Vo(T)
   par(9,1) = par(2,1) + par(4,1)*(T+298);    % KT(T)
   Kot = par(9,1);                            % for consistency in function BirchMurnsolid
   Kpt = par(3,1);
   Vo = par(8,1);
   VP(1) = BirchMurnsolid(Kot, Kpt, Vo, P);
   par(10,1) = VP(1)/1e6;                     % convert molar volume of m3/mol
   par(12,1) = par(11,1)/par(10,1);           % phase density in kg/m3

density = par(12,1);
