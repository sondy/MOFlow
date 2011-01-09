% epNd.m
% 1/4/2011
% Alessondra Springmann <asteroid@mit.edu>
% The goal: calculate epsilon 142/144 Nd for our code.

chonEp142Nd = 1.141901; % Amelin, 2004
%chonEp142Nd = 1.149498;

chonEp142NdStart = 1.1414980;

abunNd142 = 0.272; % number ratio abundance "atom percent abundance"

abunNd144 = 0.239; % number ratio abundance

abunNd142fromSm = 0.008; % number ratio abundance; see calculations in binder;
% this is from Scott (2007): ratio of 146/144 Sm = 8 * 10^-3
% [146] + [144] = 1; [146]/[144] = 8e-3
% so [146] has to be 0.008; [144] = 0.992
% all of 146 Sm has decayed to 142 Nd

massNd142 = 141.907719; % atomic mass 142-neodymium

massNd144 = 143.910083; % atomic mass 144-neodymium


% %% Just \Dpp
% 
% avgNdEER = 5.3573e-06; 
% avgSmEER = 2.4196e-06;
% 
% abunNd142 = avgNdEER*Nd142/(mass142); % from Nd
% 
% abunNd144 = avgNdEER*Nd144/(mass144); % from Nd
% 
% abunNd142fromSm = avgSmEER*Nd142fromSm/(mass142); % from Sm
% 
% ratio142_144NdSample = (abunNd142 + abunNd142fromSm)/abunNd144;
% 
% ep142Nd = ((ratio142_144NdSample/chonEp142) - 1)*10000;
% 
% fprintf('\nThe epsilon^142 Nd value for Dpp w/o residuals is %2.3g.\n\n', ...
%     ep142Nd);


%% \Dpp and residual liquids
% avgNdEERwLiq - concentration of Nd in the Dpp with residual liquids
% avgSmEERwLiq - "             "  Sm "  "   "   "    "        " 

avgNdEERwLiq = 1.8951e-05; % mass percent of Nd in the \Dpp layer
avgSmEERwLiq = 6.9755e-06; % mass percent of Sm in the \Dpp layer

% mass_res_dpp - mass of the residual liquids + the Dpp layer
% don't actually need the mass - it goes away when you divide

abunNd142wLiq = avgNdEERwLiq*abunNd142/(144.24); % from Nd

abunNd144wLiq = avgNdEERwLiq*abunNd144/(144.24); % from Nd

abunNd142fromSmwLiq = avgSmEERwLiq*abunNd142fromSm/(150.36); % from Sm

ratio142_144NdSamplewLiq = (abunNd142wLiq + abunNd142fromSmwLiq)/abunNd144wLiq;


%       1.1414980: starting 142/144 Nd ratio for earth (Carlson, 2008)
%       0.008: abundance ratio for 146Sm/144Nd (Rollinson, 2007)
ratio = 1.1414980 + 0.008*(avgNdEERwLiq/avgSmEERwLiq)*(150.36/144.24);

ep142NdwLiq = ((ratio/chonEp142Nd) - 1)*10000;

fprintf('\nThe epsilon^142 Nd value for Dpp & residuals is %2.3g.\n\n', ...
    ep142NdwLiq);

muEpDppwLiq = (ratio142_144NdSamplewLiq - chonEp142Nd)*1000000;

fprintf('\nThe mu^142 Nd value for Dpp & residuals is %2.3g.\n\n', ...
    muEpDppwLiq);

%% Emily's check: the mantle
mantleSm = 0.00001472;
mantleNd = 0.00004524;

abunNd142mantle = mantleNd*abunNd142/(massNd142);

abunNd144mantle = mantleNd*abunNd144/(massNd144);

abunNd142SmMantle = mantleSm*abunNd142fromSm/(150.36);

ratio142_144NdMantle = (abunNd142mantle + abunNd142SmMantle)/abunNd144mantle;

ep142NdMantle = ((ratio142_144NdMantle/chonEp142Nd) - 1)*10000;

fprintf('\nThe epsilon^142 Nd value for the mantle is %2.3g.\n\n', ...
    ep142NdMantle);

%% 142Nd/144Nd
disp();