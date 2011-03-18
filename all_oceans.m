% calling all oceans...

tic

globals

% depths = linspace(500, 3000, 6).*1000; 
% depths = linspace(500, 2500, 5).*1000; 
% linearly spaced magma ocean depths, all in meters

CMB_depth = 2885000;                % *** m, depth of core-mantle boundary
% from the surface of the Earth

%depths = [500000];
%depths = [2000*1000];

%%depths = [CMB_depth];

%depths = [500000, 1000000, 1500000, 2000000, 2500000];%, CMB_depth];

%%num_oceans = length(depths);

% for num_oceans = 1:length(depths)
%     ocean(depths(num_oceans))
% end

DM = CMB_depth; % depth of the magma ocean, not radius

MOFlowEarthCMB

toc