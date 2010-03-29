% calling all oceans...

global num_oceans % for making the following into vectors
global avgEER avgEERwLiq EERfracoftotal EERfracoftotalwliq avgEDR

globals

% depths = linspace(1000, 3000, 5).*1000; 
% linearly spaced magma ocean depths, all in meters

CMB = 2885000;                  % *** m, radius of core-mantle boundary
R = 6378000; 

depths = 500*1000;

%num_oceans = [];

num_oceans = 1;

five00ocean(depths)

% for num_oceans = 1:length(depths)
%     ocean(depths(num_oceans))
%end