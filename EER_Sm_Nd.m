%% EER_Sm_Nd.m
% calculates average Sm, Nd, U, Th content for EER (D") layers
% of "max" steps above CMB
% June 8, 2009

fprintf('\n \n')

max = rho2_index;
% number of steps to consider as part of the D" layer

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

%% Calculation of \Dpp REE concentrations
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

%% Residual Liquids
% these two include all the residual liquid from the surface:
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