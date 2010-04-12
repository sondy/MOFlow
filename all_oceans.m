% calling all oceans...

global num_oceans % for making the following into vectors
global avgEER avgEERwLiq EERfracoftotal EERfracoftotalwliq avgEDR

globals

% depths = linspace(500, 3000, 6).*1000; 
depths = linspace(500, 1500, 3).*1000; 
% linearly spaced magma ocean depths, all in meters

CMB = 2885000;                  % *** m, radius of core-mantle boundary
R = 6378000; 

% depths = [depths, R-CMB];

num_oceans = [];

for num_oceans = 1:length(depths)
    ocean(depths(num_oceans))
end