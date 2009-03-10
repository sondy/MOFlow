% solidusfitEarthAbe2T.m
% fits curves to solidus temperatures for Earth magma ocean
% From Abe (1997) which shows change in slope at 16 GPa

% second degree fit: S = (-0.4435)*((r/Pconversion)^2) + (64.4956)*((r/Pconversion)) + 1.2058e+003;
% third degree fit: S = (0.0282)*P^3 + (-1.9159)*P^2 + (83.5805)*P +
%       (1.1694e+003)
% fifth degree fit: S = (3.3204e-004)*((r/Pconversion)^5) + (-0.0340)*((r/Pconversion)^4) + 
%       (1.2572)*((r/Pconversion)^3) + (-20.3735)*((r/Pconversion)^2) + (182.2805)*((r/Pconversion)) + 1.1150e+003
% temperature as a function of radius:
% Tsolidus = (-1.1601e-007)*r^3 + 0.0014*r^2 + -6.3821*r + 1.4439e+004;


close all; clear all

% Don't know where these came from, but they don't quite match Abe fig. 1a
P = [130  35    30    24    20   17   16   15   10   5    0];  %GPa
T = [4200 2927  2747 2517  2337 2130 2040 2027 1907 1627 1120]; % changed T(1) from 5000


figure;
plot(P,T);
return

%R = 6378 - (P*33.3); 
R = (P - 238.5372)/-0.0374;


[Q,S] = polyfit(R,T,3);
Q
S
   
i = 0;
for k = 2900:100:6400
    i = i + 1;
    t(i) = Q(1)*k^3 + Q(2)*k^2 + Q(3)*k + Q(4);
    r(i) = k;
end

figure(1)
hold on
plot(T, R, 'x')
plot(t,r);
axis ij
xlabel('temperature');
ylabel('radius');
