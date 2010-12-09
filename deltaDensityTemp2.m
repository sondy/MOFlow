%% deltaDensityTemp2.m
% 5/5/2009 Alessondra Springmann

%% initializations

% close all;
% clear all;

% T1 = 500:100:5000; % in C
rho_2 = 3400; %kg/m^3
rho_1 = 3750; %kg/m^3

alpha = 1e-5;% 1e-5 1e-6];

%% calculation

deltaT = 500:20:2000;

max = size(deltaT,2);

deltaT = deltaT';

delta_rho = zeros(max,max);



figure(3)
hold on
for i = 1:max;
    delta_rho(i,:) = deltaT(i).*rho_1.*alpha; % change in density
    %plot(deltaT, delta_rho(i,1), '*', 'Color', [0 .5 i./max])
end


pcolor(deltaT, delta_rho(:,1), delta_rho(:,:)); %shading interp;
shading interp;
colorbar;

% for i = 1:size(T1,2);
%     delta_rho_2(i,:) = deltaT(2).*rho_1.*alpha;% - rho_1;
% end
% 
% for i = 1:size(T1,2);
%     delta_rho_3(i,:) = deltaT(3).*rho_1.*alpha;% - rho_1;
% end

%deltaT = [T2_1-T1' T2_2-T1' T2_3-T1'];

% deltaT1 = num2str(deltaT(1,1));
% deltaT2 = num2str(deltaT(1,2));
% deltaT3 = num2str(deltaT(1,3));

% deltaT1 = num2str(deltaT(1));
% deltaT2 = num2str(deltaT(2));
% deltaT3 = num2str(deltaT(3));
% 
% 
% % alpha1 = num2str(alpha(1));
% % alpha2 = num2str(alpha(2));
% % alpha3 = num2str(alpha(3));
% 
% alpha1 = num2str(alpha(1));
% alpha2 = num2str(alpha(2));
% alpha3 = num2str(alpha(3));


% legend_str = [];
% 
% legend_str1 = ['\alpha = ' alpha1 '; \DeltaT = ' deltaT1];
% legend_str2 = ['\alpha = ' alpha1 '; \DeltaT = ' deltaT2];
% legend_str3 = ['\alpha = ' alpha1 '; \DeltaT = ' deltaT3];



% plot(T1', delta_rho_2(:,1), 'r-')
% plot(T1', delta_rho_3(:,1), 'g-')
xlabel('Initial layer temperature in K')
ylabel('\Delta\rho  in kg m^{-3}')% (for \alpha = 1e-4 K^{-1})')
% legend(legend_str1, legend_str2, legend_str3, 'Location', 'East') 
hold off