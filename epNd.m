% epNd.m
% 1/4/2011
% Alessondra Springmann <asteroid@mit.edu>
% The goal: calculate epsilon 142/144 Nd for our code.

chonEp142 = 1.1414980; % Carlson, personal communication

Nd142 = 0.272;

Nd144 = 0.239;

% avgNdEERwLiq

DconNd142 = avgNdEERwLiq*Nd142;

mass142 = 141.907719;

DconNd144 = avgNdEERwLiq*Nd144;

mass144 = 143.910083;

ep142Nd = ((Nd142*mass142/Nd144*mass144)/(chonEp142 - 1))*10000;
