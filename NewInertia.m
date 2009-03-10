% Moment of inertia calculations
 
% Density of core from Bertka and Fei (1998)
% Dc = 8100 - 0.5014*r
 
% Run this after a profile, so there exists a density profile
% for the planet from the core to the surface
 
Planetmass = 6.418e23;  % mass of Mars, from NASA [kg]
 
 
% Create complete density profile for Mars
 
for m = 1:139 % For core region, assuming c-m boundary  = 1396
 
RMars(m) = (m-1)*10*1000;
DMars(m) = 8200 - (0.5)*((m-1)*10);
 
RCMars(m) = (m-1)*10*1000;
DCMars(m) = 8100 - (0.5014)*((m-1)*10);
 
end
 
for n = m + 1:m + numsteps + i
 
RMars(n) = rinv(n-m)*1000;
DMars(n) = Dsolinv(n-m);

% RMars(n) = r(n-m)*1000;
% DMars(n) = Dsol(n-m);

 
end
 
% for m = 1: k+i
% INTMars(m) = (RMars(m).^4).*DMars(m);
% end
 
Mass = 4*pi*[trapz(RMars, (RMars.^2).*DMars)];
 
Coremass = 4*pi*[trapz(RCMars, (RCMars.^2).*DCMars)]
 
I = (8/3)*pi*[trapz(RMars, (RMars.^4).*DMars)];
 
MIF = I/(Planetmass*(3396000^2))
