%% jenny craig for graphs
% 12/4/2010

%close all;

%% Cumulate density profiles
% Figures 1, 2, 3 also appear in MOFlow1 so density can be plotted before mineral recalculation that occurs in sortandinvert.m
% figure(1);title(['Density with depth for model:  ', name]); hold on; plot(D(1:marker2),r(1:marker2)/1000,'y'); xlabel('density at P and solidus temperature [kg/m3]'); ylabel('radius, km');
%figure(2);title(['Reference density with depth for model:  ', name]); hold on; plot(Do(1:marker),r(1:marker)/1000,'y'); xlabel('density at 1 atm and 1 deg C [kg/m3]'); ylabel('radius, km');
figure(3); %title(['Reference density with depth pre- and post-overturn for model:  ', name]);
hold on;
plot(Dsol, r./1000,'b.', 'LineWidth', 4)%, 'LineWidth', 4)
plot(Dsolinv, rinv./1000, 'k.', 'LineWidth', 4)%, 'LineWidth', 4); %        plot(Dsol(1:marker2), r(1:marker2)/1000,'y',...

%plot(Dsolinv(50),rinv(50)./1000,'.m', 'LineWidth', 4)   %D'' layer is
%about here

 % Now for calculated densities
density_calc_plot = zeros(maxstep, 1);

% for kk = 1:1:maxstep;
%     density_calc_plot(kk, 1) = sum(mineral_density.*liquid(kk, :)./100);
% end
% 
% plot(density_calc_plot, r./1000, 'rO');

% density_by_layer = mantle_mass_by_layer./(Mantlevolume./1000);
% % 
% % plot(density_by_layer, r./1000, 'gO');
% % plot(density_by_layer(10), r(10)./1000, 'r*');
% 
%  % Calculate density for each layer from solid composition,
%  % pre-solidification
% density_from_solid = zeros(maxstep, 1);
% 
% for kk = 1:1:maxstep;
%     density_from_solid(kk, 1) = sum(mineral_density.*solid_comp_by_layer(kk, :)./100);
% end
% 
% plot(density_from_solid, r./1000, 'cO');

xlabel('density at 1 atm and solidus temperature [kg m^{-3}]');
ylabel('radius, km');
xlim([2600 3500])
ylim([3400 6500])

legend('Location of the D" Layer', ... %...%'Reference density',...
    'Pre-overturn density', ...
    'Post-overturn density', ...
    ...%'Calculated reference density',
    ...%'Density from solid composition',...
    'Location', 'Best')

%     densityWithDepth = strcat('plots/densityWithDepth', DM_string, '.pdf');
%     print('-dpdf', densityWithDepth)
hold off;

% %% Figure 50
% figure(50);
% 
% hold on;
% 
% plot(mantle_mass_by_layer, r./1000, 'b.');
% 
% hold off;

%% liquid
% 11/15/2010

jj = 1:1:maxstep;

figure(47);

hold on;

%title('Radius versus Liquid Composition of the Coevolving Liquids');

xlabel('oxide liquid composition (mass percent)')

ylabel('radius (km)');


plot(liquid(:, 1), r./1000, 'k', 'LineWidth', 3)
plot(liquid(:, 2), r./1000, 'b', 'LineWidth', 3)
plot(liquid(:, 3), r./1000, 'c', 'LineWidth', 3)
plot(liquid(:, 4), r./1000, 'g', 'LineWidth', 3)
plot(liquid(:, 5), r./1000, 'y', 'LineWidth', 3)
% plot(liquid(:, 6), r./1000, 'r', 'LineWidth', 3)
% plot(liquid(:, 7), r./1000, 'm', 'LineWidth', 3)
% plot(liquid(:, 8), r./1000, 'Color', [0 .5 1], 'LineWidth', 3)
% plot(liquid(:, 9), r./1000, 'Color', [0 .5 .5], 'LineWidth', 3)
%plot(liquid(:, 10), r./1000, 'Color', [1 .5 0], 'LineWidth', 3)
% plot(liquid(:, 11), r./1000, 'Color', [.8 .5 1], 'LineWidth', 3)

legend('SiO_2', 'Al_2O_3', 'FeO', 'MgO', 'CaO',...
    'Location', 'EastOutside')
    %'Sm', 'Nd', 'Th','U',...
    %'OH',...% 'C',...
   % ...,
        

hold off

print('-depsc', 'plots/compositionStep.eps')

%display('foo')

% %% mantle & residual liquid mass
% % 12/4/2010
% figure(48);
% 
% hold on;
% 
% title('normalized: solidified magma ocean mass; residual liquids mass in kg')
% 
% xlabel('step number')
% 
% ylabel('normalized mass')
% 
% plot(jj', mantle_mass_vector./Mantlemass, 'k')
% plot(jj', residual_liquids_vector./Mantlemass, 'r')
% plot(jj', 0.97, 'b')
% plot(jj', 0.03, 'g')
% 
% legend('solidified magma ocean mass', 'residual liquids mass', ...
%     'Location', 'Best')
% 
% hold off;

% %% calculated density as a function of step
%
% figure(49);
%
% hold on;
%
% title('density calculated from liq_comp at each layer')
%
% xlabel('step number')
% ylabel('calculated density, kg/m^{3}')
%
% %density_calc_plot = liquid.*mineral_density;
%
%
% %
% % legend('density from liq_com'
%
% hold off;