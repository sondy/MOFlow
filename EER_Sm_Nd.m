%% EER_Sm_Nd.m
% calculates average Sm, Nd, U, Th content for EER (D") layers
% of "max" steps above CMB
% June 8, 2009

fprintf('\n \n')

max = 50; %rho2_index
% number of steps to consider as part of the D" layer

totalDvol = (4/3)*pi*(rinv(max)^3 - rinv(1)^3); % total D" volume
totalliquidvol = (4/3)*pi*(R^3 - r(maxstep)^3); 
    % total liquid unsolidified at top of MO; m^3

% preallocation
volseq = [];
delavgSm = [];
delavgNd = [];
delavgU  = [];
delavgTh = [];
    
for k = 2:max
    volseg = (4/3)*pi*(rinv(k)^3 - rinv(k-1)^3)/totalDvol;
    delavgSm(k) = solidinv(k,6)*volseg;
    delavgNd(k) = solidinv(k,7)*volseg;
    delavgU(k)  = solidinv(k,9)*volseg;
    delavgTh(k) = solidinv(k,8)*volseg;
end

%% calculation of final liquids - these might well sink with dense layers
avgNdEER = sum(delavgNd);   % for insertion into Rick Carlson's spreadsheet
avgSmEER = sum(delavgSm);
avgUEER  = sum(delavgU);
avgThEER = sum(delavgTh);

disp('For the final liquids:')
fprintf('\n')
disp(['D" Nd wt% is ', num2str(avgNdEER),...
    '; Sm wt% is ', num2str(avgSmEER)])
fprintf('\n')
disp(['U wt% is ', num2str(avgUEER),...
    '; and Th wt% is ', num2str(avgThEER), '.'])
fprintf('\n')

% write all avgEER into a global variable
avgEER(:, num_oceans) = cat(1, avgSmEER, avgNdEER, avgThEER, avgUEER);
display(avgEER)

%% these two include all the residual liquid from the surface:
% sort of like endmember of possibilities for D" composition
avgNdEERwLiq = (totalliquidvol/(totalliquidvol +...
    totalDvol))*liquid(maxstep,7) +...
    (totalDvol/(totalliquidvol+totalDvol))*avgNdEER;

avgSmEERwLiq = (totalliquidvol/(totalliquidvol +...
    totalDvol))*liquid(maxstep,6) +...
    (totalDvol/(totalliquidvol+totalDvol))*avgSmEER;

avgUEERwLiq = (totalliquidvol/(totalliquidvol +...
    totalDvol))*liquid(maxstep,9) +...
    (totalDvol/(totalliquidvol+totalDvol))*avgUEER;

avgThEERwLiq = (totalliquidvol/(totalliquidvol +...
    totalDvol))*liquid(maxstep,8) +...
    (totalDvol/(totalliquidvol+totalDvol))*avgThEER;

disp(['With all final liquids D" Nd wt% is ',num2str(avgNdEERwLiq),...
    '; Sm wt% is ', num2str(avgSmEERwLiq),...
    '; U wt% is ', num2str(avgUEERwLiq),...
    '; and Th wt% is ', num2str(avgThEERwLiq)])

fprintf('\n \n')

% write all avgEERwLiq into a global variable
avgEERwLiq(:, num_oceans) = cat(1, avgSmEERwLiq, avgNdEERwLiq,...
    avgThEERwLiq, avgUEERwLiq);

avgU = sum(delavgU);
avgTh = sum(delavgTh);
EERUfracoftotal = avgU*totalDvol/(liquid(1,9)*Mantlevolume);  
    % for comparison with Carlson's estimates
EERThfracoftotal = avgTh*totalDvol/(liquid(1,8)*Mantlevolume);    
    % of what U and Th fraction must be in D"
disp(['D" fraction of total Earth U is ', num2str(EERUfracoftotal),...
    ' and Th fraction is ', num2str(EERThfracoftotal)])
fprintf('\n')
fprintf('\n')

% write all EERfracoftotal into a global variable
EERfracoftotal(:, num_oceans) = cat(1, EERThfracoftotal, EERUfracoftotal);

%% D" fraction with liquids
EERUfracoftotalwliq = (avgU*totalDvol +...
    liquid(maxstep,9)*totalliquidvol)/(liquid(1,9)*Mantlevolume);
EERThfracoftotalwliq = (avgTh*totalDvol +...
    liquid(maxstep,8)*totalliquidvol)/(liquid(1,8)*Mantlevolume);
disp(['With all final liquids D" fraction of total Earth U is ',...
    num2str(EERUfracoftotalwliq),' and Th fraction is ',...
    num2str(EERThfracoftotalwliq)])
fprintf('\n')

% write all EERfracoftotalwliq into a global variable
EERfracoftotalwliq(:, num_oceans) = cat(1, EERThfracoftotalwliq,...
    EERUfracoftotalwliq);

avgNdEDR = mean(solidinv(max:990,7));
avgSmEDR = mean(solidinv(max:990,6));
disp(['Mantle minus D" Nd wt% is ', num2str(avgNdEDR),...
    ' and Sm wt% is ', num2str(avgSmEDR)])

fprintf('\n \n')

% write all avgNdEDR into a global variable
avgEDR(:, num_oceans) = cat(1, avgSmEDR, avgNdEDR);

%% Residual Liquid Density
densityResidual = new_mass_liquid./totalliquidvol;

%fprintf('The new mass liquid is %2.3g kg. \n', new_mass_liquid)

fprintf('The volume of the residual liquids is %2.3g kg. \n', totalliquidvol)

fprintf('The mass of the residual liquids is %2.3g kg. This is %2.3g%% the mass of the mantle. \n',...
    new_mass_liquid, (new_mass_liquid*100/mass_of_mantle))

fprintf('The density of the residual liquids is %2.3g kg/m^3 \n',...
    densityResidual)

%% Mass Balance - a test, of sorts

%fprintf('The mass of the mantle is %2.3g kg. \n', mass_of_mantle)

fprintf('The mass of the solidified ocean + residual liquids is %2.3g kg. \n',...
    (mass_solidified + new_mass_liquid))

silicateMass = (((R-DM).^3 - CMB.^3)/(R.^3 - CMB.^3))*mass_of_mantle;

% SO - solidified ocean; RL - residual liquids; US - unmelted silicates

fprintf('The mass of the SO, RL, and US is %3.4g%% of the present mantle mass. \n', ...
    ((silicateMass + mass_solidified + new_mass_liquid)*100)/(mass_of_mantle))

%% Soo, does liq_comp make sense?
% 12/4/2010

% 1: SiO2, 2: Al2O3, 3: FeO, 4: MgO, 5: CaO, 6: Sm, 7: Nd, 8: Th, 9: U, 10:
% OH, 11: C

%disp(liq_comp)

% Take the values from liq_comp and multiply them by the density of each
% species.  Does the resulting density make sense?

liq_comp_density = liq_comp.*mineral_density;

density_calc_res_liq = sum((liq_comp_density./100));

fprintf('The calculated density of the final layer is %2.3g kg/m^3. \n',...
    density_calc_res_liq)

mass_calc_res_liq = density_calc_res_liq*totalliquidvol;

% Well if it's close to 3600 kg/m^3, then yes.