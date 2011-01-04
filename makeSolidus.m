% makeSolidus.m
% March 3, 2009
% Alessondra Springmann
% Spline fit for the Abe (1997) solidus (from figure 1a in the paper)

function solidus = makeSolidus()

%% New Solidus from Abe (1997)
P = [135  35    30    24    20   17   16   15   10   5    0];  %GPa
T = [4200 2927  2747 2517  2337 2130 2040 2027 1907 1627 1120];  %K
% changed T(1) from 5000

pp_solidus = pchip(P, T); % Piecewise Cubic Hermite Interpolating Polynomial
solidus_new = @(P) ppval(pp_solidus, P);

%% Old solidus from Elkins-Tanton (2008)
pol = [-1.1601e-007, 0.0014, -6.3821, 1.4439e+004]; % from Elkins-Tanton 2008
solidus_old = @(P) polyval(pol, PtoR(P)/1000); 

%% Mix them

mix = @(P) min(P/100,1);

T = mix(P).*T + (1 - mix(P)).*solidus_old(P);
pp_solidus = pchip(P, T); % comment goes here

solidus = @(P) ppval(pp_solidus, P); %% final solidus!

% %% Test plotting
% 
%range = [0 140];

range0 = 0:1:140;

range = sort(range0, 2, 'descend');

radius_range = PtoR(range)./1000;

t_solidus_new = solidus_new(range);
t_solidus_old = solidus_old(range);
t_solidus = solidus(range);

% Plotting
% 
% figure(0); 
% 
% hold on 
% plot(t_solidus_new, radius_range, 'r', 'LineWidth', 1.5) % from Abe (1997)
% plot(t_solidus_old, radius_range, 'b', 'LineWidth', 1.5) % from Elkins-Tanton (2008)
% plot(t_solidus, radius_range, 'm', 'LineWidth', 1.5)
% legend('Abe (1997)', 'Elkins-Tanton (2008)', 'This Work',...
%     'Location', 'NorthEast')
% xlabel('Temperature (^{\circ}C)')
% ylabel('Radius (km)')
% 
% print('-depsc', 'plots/solidi.eps')


%% old comments

% at P = 0, 100% of the old solidus
% at P > 100, 100% of the new solidus

% step 1: what ratio we're going to mix these two
% step 2: mix them    
    
% P = linspace(0,140);
% mix = min(P/100,1);
% figure; plot(P,mix);  

% solidus = mix*solidus_new + (1 - mix)*solidus_old;  % HYBRID SOLIDUS

% function exists independently of its name---
% function handle to an anonymous function which takes a parameter P

% variable solidus is a handle to an anonymous function
% it doesn't have a name!
% take a parameter P and calls the matlab function ppval with pp_solidus
% (the curve fit) and P (the parameter to the anonymous function).  It then
% returns the temperature at that presure

% the beauty of this!  pp_solidus has been stored as a part of this
% function.  We can now clear pp_solidus.  As long as we keep solidus
% around, the data will still be in memory and we can access it