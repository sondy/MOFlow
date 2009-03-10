% graphsdeep

disp('In graphsdeep')

%% Cumulate density profiles
% Figures 1, 2, 3 also appear in MOFlow1 so density can be plotted before mineral recalculation that occurs in sortandinvert.m
% figure(1);title(['Density with depth for model:  ', name]); hold on; plot(D(1:marker2),r(1:marker2)/1000,'y'); xlabel('density at P and solidus temperature [kg/m3]'); ylabel('radius, km');
%figure(2);title(['Reference density with depth for model:  ', name]); hold on; plot(Do(1:marker),r(1:marker)/1000,'y'); xlabel('density at 1 atm and 1 deg C [kg/m3]'); ylabel('radius, km');
figure(3); title(['Reference density with depth pre- and post-overturn for model:  ', name]); 
    hold on; 
    plot(Dsol(1:marker2), r(1:marker2)/1000,'y',...
        Dsolinv, rinv/1000, 'k'); 
    xlabel('density at 1 atm and solidus temperature [kg m^{-3}]'); 
    ylabel('radius, km');    
%    print -depsc 'plots/densityWithDepth.eps'
    
% figure(4); title(['Phase densities for model:  ', name]); hold on; axis([2000 4500 5400 6400]); plot(perov, r/1000, magnesio, r/1000, bet, r/1000, gam, r/1000, maj, r/1000, stish, r/1000, cpyrox, r/1000, opyrox, r/1000, gar, r/1000, spi, r/1000, plag, r/1000, oliv, r/1000, iliq, r/1000); 
%     legend('perovskite', 'magnesiowustite', 'beta', 'gamma', 'majorite', 'stishovite', 'cpx', 'opx', 'garnet', 'spinel', 'plagioclase', 'olivine', 'int liq', 0); xlabel('Density [kg/m3]'); ylabel('radius, km'); zoom on;
% figure(5); title(['Phase densities for model:  ', name]); axis ij; hold on; plot(perov, P, magnesio, P, bet, P, gam, P, maj, P, stish, P, cpyrox, P, opyrox, P, gar, P, spi, P, plag, P, oliv, P, iliq, P); 
%     legend('perovskite', 'magnesiowustite', 'beta', 'gamma', 'majorite', 'stishovite', 'cpx', 'opx', 'garnet', 'spinel', 'plagioclase', 'olivine', 'int liq', 0); xlabel('Density [kg/m3]'); ylabel('Pressure [GPa]'); zoom on; %axis([2000 4500 5400 6400]);
%figure(5); title(['Solidus density with depth for overturned model']);  hold on;
%    plot(Dsolinv,rinv/1000, 'r',Dinv,rinv/1000, 'y', Doinv,rinv/1000, 'k'); xlabel('density [kg/m3]'); ylabel('radius, km'); legend('Density at 1 atm, solidus','Density at P, solidus','Density at 1 atm, 1 deg',0);
 

%% Temperatures and pressures
figure(6); title(['temperature evolutions [C] for: ', name]); hold on; plot(time/3.1536e13, Tsurf, 'r-', time/3.1536e13, Tsolid, 'b-'); legend('surface temp','interior temp',0); xlabel('time [Ma]')
%%    figure(7); title(['temperature evolutions [C] for: ', name]); hold on; plot(Tsurffig, R/1000, 'ro', Tsolidfig, rfig/1000, 'bx', Tsolid(1) - adiabslope*(rfig - (R-RM)), rfig/1000, 'k-'); legend('surface temp','temp at top of solid','adiabat',0); ylabel('radius [km]')
figure(8); title(['Temperature with depth for overturned model']);  hold on;
    plot(Tinv,rinv/1000, 'r', CorrectedTinv,rinv/1000, 'g', Tsolid, r/1000, 'k:'); xlabel('temperature [C]'); ylabel('radius, km'); legend('Inverted temperature','Adiabatically corrected','solidus',0);
%    print -depsc 'plots/temperatureWithDepth.eps'
%   plot(Tinv,rinv/1000, 'r', CorrectedTinv,rinv/1000, 'g', solidus, rinv/1000, 'k:'); xlabel('temperature [C]'); ylabel('radius, km'); legend('Inverted temperature','Adiabatically corrected','solidus',0);
%figure(9); title(['Pressure vs r for model:  ', name]); hold on; plot(P,r/1000, 'r'); xlabel('Pressure [GPa]'); ylabel('radius, km');
%%    figure(9); title(['temperature evolutions [C] for: ', name]); hold on; plot(Tsurf, Patm*10^-5, 'r-'); legend('surface temp',0); xlabel('Surface T [C]'); ylabel('Patm [bar]')
%%    figure(10); title(['temperature evolutions [C] for: ', name]); hold on; plot(Tsurf, log10(HPatm*10^-5), 'r-'); legend('surface temp',0); xlabel('Surface T [C]'); ylabel('log10(HPatm) [bar]')

% Flux, emissivity, solidification rate, viscous boundary layer
%figure(10); title(['log(10)planetary heat flux [W/m^2] for: ', name]); hold on; plot(time/3.1536e13, log10(flux), 'g-'); xlabel('time [Ma]'); %axis([0 1.2 0.5 4])
%figure(11); title(['log(10)emissivity for: ', name]); hold on; plot(time/3.1536e13, log10(emiss), 'm-'); xlabel('time [Ma]')
%figure(12); title(['Convective velocity [m/year] for: ', name]); hold on; plot(time/3.1536e13, veloc*3.1536e7, 'r-'); xlabel('time [Ma]')
%figure(12); title(['Atmospheric water partial pressure vs tau*: ', name]); hold on; plot(HPatm, taustarw, 'k-'); xlabel('Pressure [Pa]'); ylabel('Tau star')
%figure(13); title(['Viscous boundary layer [m] for: ', name]); hold on; plot(time/3.1536e13, delt, 'r-'); xlabel('time [Ma]')
%figure(13); title(['Flux [W/m^2] vs surface temperature: ', name]); hold on; plot(Tsurf, log10(flux), 'k-'); xlabel('T [C]'); ylabel('Log(10)Flux [W/m^2]')
figure(14); title(['Solidification rate by normalized volume: ', name]); hold on; plot(time/3.1536e13, vol/1000, 'k-'); xlabel('time [Ma]')

%% Atmospheres
% figure(15); title(['log(10)mass added to atmosphere [kg] for: ', name]); hold on; plot(time/3.1536e13, log10(Catmadd), 'ro', time/3.1536e13, log10(Hatmadd), 'bx'); legend('carbon dioxide','water',0); xlabel('time [Ma]')
% figure(16); title(['Atmospheric pressure and partial pressures [bar] for: ', name]); hold on; plot(time/3.1536e13, Patm*10^-5, 'k-', time/3.1536e13, CPatm*10^-5, 'r-', time/3.1536e13, HPatm*10^-5, 'b-'); legend('total', 'carbon species', 'water',0); xlabel('time [Ma]')
% figure(17); title(['Atmospheric pressure [bar] for: ', name]); hold on; plot(Tsurf, CPatm*10^-5, 'r-', Tsurf, HPatm*10^-5, 'b-'); legend('carbon dioxide','water',0); xlabel('surface temp [C]'); %axis([0 1.2 0 3e21])
% figure(18); title(['fractions of total volatiles into atmosphere for model: ', name]); hold on; plot(time/3.1536e13, Cfractionout, 'r-', time/3.1536e13, Hfractionout, 'b-'); legend('carbon degassing','water degassing'); xlabel('time [Ma]')

%% Magma ocean liquids
%figure(19); title(['Liquid water content [mass percent] for: ', name]); hold on; plot(time/3.1536e13, H2Osaturation, 'r-', time/3.1536e13, liquid(:,10), 'b-'); legend('saturation','magma content',0); xlabel('time [Ma]')
%%    figure(20); title(['Liquid water content [mass percent] for: ', name]); hold on; plot(Tsurf, liquid(:,10)); xlabel('surface temperature')
%figure(20); title(['Liquid carbon content [mass percent] for: ', name]); hold on; plot(time/3.1536e13, CO2saturation, 'r-', time/3.1536e13, liquid(:,11), 'b-'); legend('saturation','magma content',0); xlabel('time [Ma]'); %axis([0 0.9 0 3])
%figure(21); title(['Volatile contents in evolving liquid']); hold on; plot(liquid(:,10),r/1000, 'g', liquid(:,11),r/1000,'r'); xlabel('water and carbon content in evolving liquid wt%]'); ylabel('radius, km'); legend('water','carbon',0);
%figure(22); title(['Oxide percentage of evolving liquid:  ', name]); hold on; plot(liquid(:,1), r/1000, liquid(:,2), r/1000, liquid(:,3), r/1000, liquid(:,4), r/1000, liquid(:,5), r/1000); xlabel('wt% of oxide'); ylabel('radius, km'); legend('SiO2','Al2O3', 'Feo','MgO', 'CaO', 0);
% figure(23); title(['Ppm of trace elements in evolving liquids for model: ', name]); hold on;
%    plot(liquid(:,6)*10^4, r/1000, liquid(:,7)*10^4, r/1000, liquid(:,8)*10^4, r/1000, liquid(:,9)*10^4, r/1000); xlabel('ppm of trace element'); ylabel('radius, km'); legend('Sm','Nd','Lu','Hf',0);

%% Major element compositions of cumulates
%%    figure(24); title(['Oxide percentage with depth for pre-overturn model:  ', name]); hold on; plot(solid(:,1), r/1000, solid(:,2), r/1000, solid(:,3), r/1000, solid(:,4), r/1000, solid(:,5), r/1000); xlabel('wt% of oxide'); ylabel('radius, km'); legend('SiO2','Al2O3', 'Feo','MgO', 'CaO', 0);
%figure(25); title(['Oxide percentage with depth for overturned model:  ', name]); hold on;
%   plot(solidinv(:,1), rinv/1000, solidinv(:,2), rinv/1000, solidinv(:,3), rinv/1000, solidinv(:,4), rinv/1000, solidinv(:,5), rinv/1000); xlabel('wt% of oxide'); ylabel('radius, km'); legend('SiO2','Al2O3','FeO','MgO','CaO',0);
%%    figure(26); title(['Mg# of phases for pre-overturn model:  ', name]); hold on; plot(MG(:,1),r/1000, MG(:,2),r/1000); xlabel('Mg#'); ylabel('radius, km');legend('Liquid','Bulk solid',0);
% figure(27); title(['Mg# of overturned cumulates']); hold on; plot(MGinv(:,2),rinv/1000, 'g'); hold on; xlabel('Mg#'); ylabel('radius, km'); axis([0 1 1396 3396]);

%% Volatile content of cumulates
%figure(28); title(['Mass of volatiles in cumulates: ', name]); hold on; plot(time/3.1536e13, solid(:,11), 'r-', time/3.1536e13, solid(:,10), 'b-', time/3.1536e13, solid(:,10)+ solid(:,11), 'k-'); legend('carbon dioxide','water','total',0); xlabel('time [Ma]')
%figure(29); title(['Water content in cumulates for model: ', name]); hold on; plot(solidinv(:,10),rinv/1000, 'r', solid(:,10),rinv/1000, 'k'); xlabel('water content in cumulates wt%]'); ylabel('radius, km'); legend('after overturn','before overturn',0);
%figure(30); title(['Carbon content in cumulates for model: ', name]); hold on; plot(solidinv(:,11),rinv/1000, 'r', solid(:,11),rinv/1000, 'k'); xlabel('carbon content in cumulates wt%]'); ylabel('radius, km'); legend('after overturn','before overturn',0);
%figure(30); title(['Water content in cumulates for model: ', name]); hold on; plot(solidOHinv(:),rinv/1000, 'r', intliqOHinv(:),rinv/1000, 'b', solidOHinv(:)+intliqOHinv(:),rinv/1000, 'k'); xlabel('water content in cumulates wt%]'); ylabel('radius, km'); legend('in nominally anhydrous minerals','in interstitial liquids','total',0);

%%
% Trace element content of cumulates
% figure(31); title(['ppm elements in solids for model:  ', name]); hold on; 
%    plot(solid(:,6)*10^4, r/1000, solid(:,7)*10^4, r/1000, solid(:,8)*10^4/0.0243, r/1000, solid(:,9)*10^4/0.1040, r/1000); xlabel('ppm'); ylabel('radius, km');legend('Sm','Nd','Lu','Hf',0);
% figure(32); title(['Sm/Nd trace element ratios CI chondrite normalized for pre-overturn model:  ', name]); axis([1 8 5800 6400]);
%      hold on; plot((solid(:,6)*10^4/0.1471)./(solid(:,7)*10^4/0.4524),r/1000,'r-'); %,(solid(:,8)*10^4/0.0243)./(solid(:,9)*10^4/0.1040), r/1000, 'k-'); 
%      xlabel('chondrite-normalized ratios'); ylabel('radius, km');
%      legend('Source (Sm/Nd)N','Source (Th/U)N');
% figure(33); title(['Trace element ratios CI chondrite normalized for overturned model:  ', name]);
%      hold on; plot((solidinv(:,6)*10^4/0.1471)./(solidinv(:,7)*10^4/0.4524),...
%          (solidinv(:,8)*10^4/0.0243)./(solidinv(:,9)*10^4/0.1040), 'rs',...
%          (liquid(:,6)*10^4/0.1471)./(liquid(:,7)*10^4/0.4524),...
%          (liquid(:,8)*10^4/0.0243)./(liquid(:,9)*10^4/0.1040), 'go');
%      xlabel('Source (Sm/Nd)N'); 
%      ylabel('Source (Th/U)N'); 
     %axis([0.8 1.5 0.5 3.0]);
figure(34); title(['Trace element concentrations chondrite-normalized with depth for overturned model:  ', name]);
     hold on; plot(solidinv(:,6)*10^4/0.1471, rinv/1000, solidinv(:,7)*10^4/0.4524,... 
         rinv/1000, solidinv(:,8)*10^4/0.0243, rinv/1000, solidinv(:,9)*10^4/0.1040, rinv/1000); 
     xlabel('fraction of chondritic'); 
     ylabel('radius, km'); 
     legend('Sm','Nd','Th','U',0); %axis([0 .5 2000 3400]);
%figure(35); title(['Sm/Nd reservoir ratio with depth for overturned model:  ', name]); hold on; plot(SmNdinv,rinv/1000); xlabel('147Sm/144Nd initial ratio of reservoir'); ylabel('radius [km]');
% figure(36); title(['176Lu/177Hf for overturned model:  ', name]); hold on; plot(Lu6Hf7inv, rinv/1000); xlabel('176Lu/177Hf'); ylabel('radius, km');

%% Crust formation
%figure(37); title(['Total melt ', num2str(totalmelt/(10^9)),' km^3, covers planet to ',num2str(depth),' km']);
%    hold on; plot(meltpercenttop,rinv/1000, meltpercentbottom, rinv/1000); xlabel('percentage melt created by overturn'); ylabel('radius [km]'); axis([0 100 1396 3396]);

% FIGURE 38 resides in script "sortandinvert.m"

%% Cool to clement
%figure(39); title(['model ', name,' for ', num2str(tfinal/(3.14*10^13)),' Ma']); plot(temp(:,1)-273, r2c, 'k-', halftemp(:)-273, r2c, 'm-', newtemp(:)-273, r2c,'g-'); legend('initial temperature profile','midway temperature profile', 'final temperature profile'); xlabel('temperature [C]'); ylabel('radius')
%figure(39); title(['Temp profiles for ', name,' for ', num2str(tfinal/(3.14*10^13)),' Ma']); plot(temp(:,1)-273, r2c, 'k-', quartertemp-273, r2c, 'r-', halftemp-273, r2c, 'r-', threeqtemp-273, r2c, 'r-', tempfinal(:)-273, r2c,'g-'); legend('initial temperature profile', '1/4 time', '1/2 time', '3/4 time', 'final temperature profile',0); xlabel('temperature [C]'); ylabel('radius')
    %axis([800 3000 1396 3396]);
%figure(40); title(['Surface heat flux for ', name,' for ', num2str(tfinal/(3.14*10^13)),'Ma']); plot(time2c, surfaceflux, 'r-'); xlabel('time [Ma]'); ylabel('flux [J/m2sec]');
%    figure(41); title(['Surface temp ', name,' for ', num2str(tfinal/(3.14*10^13)),'Ma']); plot(time2c, temp(992,:), 'g-', time2c, temp(991,:), 'r:'); legend('top b.c. temp','surface temp',0); xlabel('time [Ma]'); ylabel('temp [C]');


% x1 = [0:.1:40];
% y1 = 4.*cos(x1)./(x1+2);
% x2 = [1:.2:20];
% y2 = x2.^2./x2.^3;
% hl1 = line(x1,y1,'Color','r');
% ax1 = gca;
% ax2 = axes('Position',get(ax1,'Position'),...
%            'XAxisLocation','top',...
%            'YAxisLocation','right',...
%            'Color','none',...
%            'XColor','k','YColor','k');
% hl2 = line(x2,y2,'Color','k','Parent',ax2);
% 
% hold on;

figure(42);

hold on;

x1 = Dsolinv;
y1 = rinv/1000;

x2 = solidinv(:,6)*10^4/0.1471;
y2 = rinv/1000;
x3 = solidinv(:,7)*10^4/0.4524;
y3 = rinv/1000; 
x4 = solidinv(:,8)*10^4/0.0243; 
y4 = rinv/1000; 
x5 = solidinv(:,9)*10^4/0.1040; 
y5 = rinv/1000;

hl1 = line(x1, y1, 'Color', 'k');
xlabel('density at 1 atm and solidus temperature [kg/m^3]');
ylabel('radius, km');
ax1 = gca;

ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
hl2 = line(x2,y1,'Color','b','Parent',ax2);
hl3 = line(x3,y1,'Color','g','Parent',ax2);
hl4 = line(x4,y1,'Color','r','Parent',ax2);
hl5 = line(x5,y1,'Color','c','Parent',ax2);

legend('Sm','Nd','Th','U',0); %axis([0 .5 2000 3400]);

title(['Trace element concentrations chondrite-normalized with depth for overturned model:  ', name]);

hold off;


% 
% figure(42);...
% %    title(['Reference density with depth pre- and post-overturn for model:  ', name]);...
%     plot(...%Dsol(1:marker2),r(1:marker2)/1000,'y',
%     Dsolinv,rinv/1000, 'k');...
%     xlabel('density at 1 atm and solidus temperature [kg/m3]');...
%     ylabel('radius, km');
% 
% 
% 
% %figure(43); 
% %title(['Trace element concentrations chondrite-normalized with depth for overturned model:  ', name]); 
% plot(solidinv(:,6)*10^4/0.1471, rinv/1000, solidinv(:,7)*10^4/0.4524, rinv/1000, solidinv(:,8)*10^4/0.0243, rinv/1000, solidinv(:,9)*10^4/0.1040, rinv/1000);...
%     xlabel('fraction of chondritic');...
%     ylabel('radius, km');...
%     legend('Sm','Nd','Lu','Hf',0); %axis([0 .5 2000 3400]);

