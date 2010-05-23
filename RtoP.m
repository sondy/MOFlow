% RtoP.m
% March 3, 2009
% Alessondra Springmann
% Takes the distance from the center of the Earth in meters and converts 
% it into pressure in GPa

function pressure = RtoP(radius)

% linear relationship from MOFlow code
CMB_depth = 2885000;             % *** m, depth to core-mantle boundary
R = 6378000;                     % *** m, total radius of planet
CMB = R - CMB_depth;

m = 140*1000/(CMB - R);

b = 140 - m*(CMB/1000);

pressure = m*(radius/1000) + b; 

