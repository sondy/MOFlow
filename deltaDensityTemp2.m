%% deltaDensityTemp2.m
% 5/5/2009 Alessondra Springmann

%% initializations

% close all;
% clear all;

% T1 = 500:100:5000; % in C
rho_2 = 3400; %kg/m^3
rho_1 = 3300; %kg/m^3

alpha = 1e-5;% 1e-5 1e-6];

%% calculation

deltaT = 500:50:5000;

max = size(deltaT,2);

deltaT = deltaT';

delta_rho = zeros(max,max);

%figure(3)
hold on
for i = 1:max;
    delta_rho(i,:) = deltaT(i).*rho_1.*alpha; % change in density
    %plot(deltaT, delta_rho(i,1), '*', 'Color', [0 .5 i./max])
end

pcolor(deltaT, delta_rho(:,1), delta_rho(:,:));
shading interp;
colorbar;
colormap('hot')

xlabel('\Delta T in K')
ylabel('\Delta\rho  in kg\cdotm^{-3} (for \alpha = 1e-5 K^{-1})')
% legend(legend_str1, legend_str2, legend_str3, 'Location', 'East') 
hold off

print('-dpng', 'plots/deltaRho.png')