% globals.m
% June 7, 2009 Alessondra Springmann
% a bunch of initializations to make running MOFlowEarth* be more pleasant

%global num_oceans % for making the following into vectors
%global avgEER avgEERwLiq EERfracoftotal EERfracoftotalwliq avgEDR
%global all_liquid_composition % keep tabs on what the liquid composition
%is doing!

global Dolinv liquid maxstep num_oceans r rinv solidinv MW
% for later plotting/debugging
%global num_oceans % for making the following into vectors
global densityResidual
global avgEER avgEERwLiq EERfracoftotal EERfracoftotalwliq avgEDR
%global KSm_m KNd_m KTh_m KU_m KOH_m KCO_m


% global xrecord xrecordn % for determing ranges for fzero in BMsolid.m
% 
% xrecordn = 0; % when you initialize a global, it gets cast as a vector.  
% % Incrementing a vector results in an empty vector.  Who knew.