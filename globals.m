% globals.m
% June 7, 2009 Alessondra Springmann
% a bunch of initializations to make running MOFlowEarth* be more pleasant

global Dolinv liquid maxstep r rinv solidinv % for later plotting/debugging
global xrecord xrecordn % for determing ranges for fzero in BMsolid.m
global num_oceans % for making the following into vectors
global avgEER avgEERwLiq EERfracoftotal EERfracoftotalwliq avgEDR

xrecordn = 0; % when you initialize a global, it gets cast as a vector.  
% Incrementing a vector results in an empty vector.  Who knew.  Thanks
% Henry.

