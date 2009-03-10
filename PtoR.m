% PtoR.m
% March 3, 2009
% Alessondra Springmann
% takes a pressure in GPa and converts it to a distance from the center of 
% the Earth in m

function radius = PtoR(pressure)

radius = (1000/-0.0374)*(pressure -238.5372);

% p =  -0.0374*(radius/1000) + 238.5372; 
% 
% 1000*(p - 238.5372)/-0.0374