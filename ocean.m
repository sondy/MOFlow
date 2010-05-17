%clear all;

function ocean(depth)

globals

DM = depth; % depth of the magma ocean, not radius

tic

% global xrecord xrecordn
% xrecord=zeros(1,100000);
% xrecordn=0;

MOFlowEarthCMB

toc