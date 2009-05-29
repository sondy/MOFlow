totalliquidvol;

% 1: SiO2, 2: Al2O3, 3: FeO, 4: MgO, 5: CaO
final_liq_perc = [2.5369, 20.2142, 26.634, 0.9948, 49.6088];

final_liq_rho = [2.24, 2.70, 5.75, 3.62, 4.5]; % from Grove & Baker 1983

rho_weighted = final_liq_perc.*final_liq_rho./100;

rho_ave = sum(rho_weighted)