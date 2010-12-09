%% jenny craig for graphs
% 12/4/2010

%close all;

%% Cumulate density profiles
% Figures 1, 2, 3 also appear in MOFlow1 so density can be plotted before mineral recalculation that occurs in sortandinvert.m
% figure(1);title(['Density with depth for model:  ', name]); hold on; plot(D(1:marker2),r(1:marker2)/1000,'y'); xlabel('density at P and solidus temperature [kg/m3]'); ylabel('radius, km');
%figure(2);title(['Reference density with depth for model:  ', name]); hold on; plot(Do(1:marker),r(1:marker)/1000,'y'); xlabel('density at 1 atm and 1 deg C [kg/m3]'); ylabel('radius, km');
figure(3); %title(['Reference density with depth pre- and post-overturn for model:  ', name]);
hold on;
plot(Dsol, r./1000,'bO')%, 'LineWidth', 4)
plot(Dsolinv, rinv./1000, 'kO')%, 'LineWidth', 4); %        plot(Dsol(1:marker2), r(1:marker2)/1000,'y',...

plot(Dsolinv(50),rinv(50)./1000,'Om')   %D'' layer is about here

 % Now for calculated densities
density_calc_plot = zeros(maxstep, 1);

% for kk = 1:1:maxstep;
%     density_calc_plot(kk, 1) = sum(mineral_density.*liquid(kk, :)./100);
% end
% 
% plot(density_calc_plot, r./1000, 'rO');

density_by_layer = mantle_mass_by_layer./(Mantlevolume./1000);
% 
% plot(density_by_layer, r./1000, 'gO');
% plot(density_by_layer(10), r(10)./1000, 'r*');

 % Calculate density for each layer from solid composition,
 % pre-solidification
density_from_solid = zeros(maxstep, 1);

for kk = 1:1:maxstep;
    density_from_solid(kk, 1) = sum(mineral_density.*solid_comp_by_layer(kk, :)./100);
end

plot(density_from_solid, r./1000, 'cO');

xlabel('density at 1 atm and solidus temperature [kg m^{-3}]');
ylabel('radius, km');
xlim([2600 3500])
ylim([3400 6500])

legend(...%'Reference density',...
    'Pre-overturn density',...
    'Post-overturn density', 'Location of the D" Layer', ...
    ...%'Calculated reference density',
    'Density from solid composition',...
    'Location', 'Best')

%     densityWithDepth = strcat('plots/densityWithDepth', DM_string, '.pdf');
%     print('-dpdf', densityWithDepth)
hold off;

%% Figure 50
figure(50);

hold on;

plot(mantle_mass_by_layer, r./1000, 'b.');

hold off;

%% liquid
% 11/15/2010

jj = 1:1:maxstep;

figure(47);

hold on;

title('Liquid Composition as a function of step');

xlabel('step number')

ylabel('liquid composition');

plot(jj', liquid(:, 1), 'k');
plot(jj', liquid(:, 2), 'b');
plot(jj', liquid(:, 3), 'c')
plot(jj', liquid(:, 4), 'g')
plot(jj', liquid(:, 5), 'y')
plot(jj', liquid(:, 6), 'r')
plot(jj', liquid(:, 7), 'm')
plot(jj', liquid(:, 8))
plot(jj', liquid(:, 9))
plot(jj', liquid(:, 10))
plot(jj', liquid(:, 11))

legend('SiO2', 'Al2O3', 'FeO', 'MgO', 'CaO', 'Sm', 'Nd', 'Th', 'U',...
    'OH', 'C', 'Location', 'Best')

hold off

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