% RtoP.m
% March 3, 2009
% Alessondra Springmann
% Takes the distance from the center of the Earth in meters and converts 
% it into pressure in GPa

function pressure = RtoP(radius)

% linear relationship from MOFlow code

pressure = -0.0374*(radius/1000) + 238.5372; 

