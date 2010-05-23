% PtoR.m
% March 3, 2009
% Alessondra Springmann
% takes a pressure in GPa and converts it to a distance from the center of 
% the Earth in m

function radius = PtoR(pressure)

CMB_depth = 2885000;             % *** m, depth to core-mantle boundary
R = 6378000;                     % *** m, total radius of planet
CMB = R - CMB_depth;

m = 140*1000/(CMB - R);

b = 140 - m*(CMB/1000);

radius = (1000/m)*(pressure - b);

% p =  -0.0374*(radius/1000) + 238.5372; 
% 
% 1000*(p - 238.5372)/-0.0374