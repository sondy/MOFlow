% calling all oceans...

global num_oceans % for making the following into vectors
global avgEER avgEERwLiq EERfracoftotal EERfracoftotalwliq avgEDR
%global all_liquid_composition % keep tabs on what the liquid composition is doing!

globals

% depths = linspace(500, 3000, 6).*1000; 
depths = linspace(500, 2500, 5).*1000; 
% linearly spaced magma ocean depths, all in meters

CMB_depth = 2885000;                % *** m, depth of core-mantle boundary
% from the surface of the Earth

%depths = [CMB_depth];

%depths = [500000];

%depths = 2000*1000;

%depths = [depths]
%depths = [depths, CMB_depth];
%depths = [2000*1000];
%depths = [CMB_depth];

num_oceans = [];

%all_liquid_composition = [];

for num_oceans = 1:length(depths)
    ocean(depths(num_oceans))
end