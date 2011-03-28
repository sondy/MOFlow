% plotMonte
% 3/27/2011
% Alessondra Springmann

% takes a .dat file and plots interesting things

%% initializations

load runInfo.dat;

Ca_Kd_Sm = runInfo(:,8);
Ca_Kd_Nd = runInfo(:,9);

MgFe_Kd_Sm = runInfo(:,14);
MgFe_Kd_Nd = runInfo(:,15);

EDR_mu = runInfo(:,20);
EER_mu = runInfo(:,21);

%% Ca_Sm
figure(90);
hold on
scatter(EDR_mu, EER_mu, 5, Ca_Kd_Sm)
plot(0, -54, 'ks')
colorbar
colormap(copper)
xlabel('EDR \mu^{142}Nd (ppm)')
ylabel('EER \mu^{142}Nd (ppm)')
hold off

print -depsc plots/runInfo_Ca_Sm.eps

%% Ca_Nd
figure(91);
hold on
scatter(EDR_mu, EER_mu, 5, Ca_Kd_Nd)
plot(0, -54, 'ks')
colorbar
colormap(copper)
xlabel('EDR \mu^{142}Nd (ppm)')
ylabel('EER \mu^{142}Nd (ppm)')
hold off

print -depsc plots/runInfo_Ca_Nd.eps

%% MgFe_Sm
figure(92);
hold on
scatter(EDR_mu, EER_mu, 5, MgFe_Kd_Sm)
plot(0, -54, 'ks')
colorbar
colormap(copper)
xlabel('EDR \mu^{142}Nd (ppm)')
ylabel('EER \mu^{142}Nd (ppm)')
hold off

print -depsc plots/runInfo_MgFe_Sm.eps

%% MgFe_Nd
figure(93);
hold on
scatter(EDR_mu, EER_mu, 5, MgFe_Kd_Nd)
plot(0, -54, 'ks')
colorbar
colormap(copper)
xlabel('EDR \mu^{142}Nd (ppm)')
ylabel('EER \mu^{142}Nd (ppm)')
hold off

print -depsc plots/runInfo_Ca_Sm.eps