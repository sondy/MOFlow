% makeSolidus.m
% September 23, 2009, Sonia Tikoo

function solidus = soniaSolidus()

% **solidus from Abe (1997)** 
%   this is a set of individual P,T values:

P = [130  35    30    24    20   17   16   15   10   5    0];  %GPa
Tabe = [4200 2927  2747 2517  2337 2130 2040 2027 1907 1627 1120]; % changed T(1) from 5000

figure(100);
plot(P, Tabe, 'r');

hold on 

R = (P - 238.5372)/-0.0374;
range = [0 140];
Pplot = 0:0.1:140;
Rplot = (Pplot - 238.5372)/-0.0374;
% **solidus from Elkins-Tanton (2008)** 
%   this is a 3rd degree polynomial calculating T from a given R:
%   Tearth = (-1.160e-7)*R^3 + (0.0014)*R^2 - (6.382)*R + 1.444e4

% This line uses Lindy's funtion to evaluate T at Abe's pressures

Telk = (-1.160e-7)*R.^3 + (0.0014)*R.^2 - (6.382)*R + 1.444e4;
TelkPlot = (-1.160e-7)*Rplot.^3 + (0.0014)*Rplot.^2 - (6.382)*Rplot + 1.444e4;
plot(Pplot, TelkPlot, 'b');

Tmin = min(Tabe,Telk);

% plot(P, Tmin, 'g')

pp_solidus_sonia = spline(P,Tmin); % returns piecewise polynomial form of cubic spline
% pp_solidus_elkins = spline(P,Telk); % returns piecewise polynomial form of cubic spline
% pp_solidus_abe = spline(P,Tabe); % returns piecewise polynomial form of cubic spline

% print out the coeficients

pp_solidus_sonia.coefs;


% range = [0 140];
% 
% figure
% hold on 
% fplot(pp_solidus_abe, range, 'r') % from Abe (1997)
% fplot(pp_solidus_elkins, range, 'b') % from Elkins-Tanton (2008)
% fplot(pp_solidus_sonia, range, 'g')
% legend('New', 'Old', 'Hybrid', 'Location', 'NorthWest')
% xlabel('Pressure in GPa')
% ylabel('Temperature in C')
% 


solidus = @(P) ppval(pp_solidus_sonia, P); %% final solidus!

fplot(solidus, range, 'g');


return





