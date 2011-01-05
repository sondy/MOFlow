% epNd.m
% 1/4/2011
% Alessondra Springmann <asteroid@mit.edu>
% The goal: calculate epsilon 142/144 Nd for our code.

chonEp142 = 1.1414980; % Carlson, personal communication

Nd142 = 0.272;

Nd144 = 0.239;

Nd142fromSm = 0.008; % see calculations in binder;
% this is from Scott (2007): ratio of 146/144 Sm = 8 * 10^-3
% [146] + [144] = 1; [146]/[144] = 8e-3
% so [146] has to be 0.008; [144] = 0.992
% all of 146 Sm has decayed to 142 Nd

mass142 = 141.907719;

mass144 = 143.910083;

% avgNdEERwLiq - concentration of Nd in the Dpp with residual liquids
% avgSmEERwLiq - "             "  Sm "  "   "   "    "        " 

avgNdEERwLiq = 1.8951e-05;
avgSmEERwLiq = 6.9755e-06;

% mass_res_dpp - mass of the residual liquids + the Dpp layer
% don't actually need the mass - it goes away when you divide

abunNd142 = avgNdEERwLiq*Nd142/(mass142); % from Nd

abunNd144 = avgNdEERwLiq*Nd144/(mass144); % from Nd

abunNd142fromSm = avgSmEERwLiq*Nd142fromSm/(mass142); % from Sm

ratio142_144NdSample = (abunNd142 + abunNd142fromSm)/abunNd144;

ep142Nd = ((ratio142_144NdSample/chonEp142) - 1)*10000;

fprintf('The epsilon^142 Nd value for our model is %2.3g.  \n', ...
    ep142Nd);