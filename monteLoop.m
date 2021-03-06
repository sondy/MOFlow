% monteLoop.m
% Alessondra Springmann, 3/19/2011
% Monte Carlo simulations for MOFlow, specifically, changing the partition
% coefficients for perovskite, post-perovskite, and magnesiowustite

% decide how many runs will occur
% create a range of variables for each run
% write run number, input variables into a file
% calculate output values, write into a file

%matlabpool local 4

tic

numRuns = 1;

numPar = 6;
numMin = 3;
kdNum = numPar*numMin; % number of partition coefficients * number of minerals

% %% magnesiowustite
% magSmNdUMin = 0.0001;
% magSmNdUMax = 0.1;
% 
% magThMin = 0.001;
% magThMax = 1.0;
% 
% magOHMin = 0.0001;
% magOHMax = 0.1;
% 
% magCOMin = 0.00001;
% magCOMax = 0.01;
% 
% Kd_m_minmax = [magSmNdUMin, magSmNdUMax;...
%     magSmNdUMin, magSmNdUMax;...
%     magThMin, magThMax ;...
%     magSmNdUMin, magSmNdUMax;...
%     magOHMin, magOHMax;...
%     magCOMin, magCOMax];

%% perovskite, Ca

% 1 - 20 for Sm
% 5 - 20 for Nd
% 0.01 - 20 for U
% 1 - 20 for Th

perCaSmMin = 1;
perCaSmMax = 5;

perCaNdMin = 1;
perCaNdMax = 10;

perCaUMin = 0.01;
perCaUMax = 1;

perCaThMin = 1;
perCaThMax = 10;

perCaOHMin = 0.00001;
perCaOHMax = 0.001;

perCaCOMin = 0.00005;
perCaCOMax = 0.005;

Kd_p_Ca_minmax = [perCaSmMin, perCaSmMax;...
    perCaNdMin, perCaNdMax;...
    perCaThMin, perCaThMax ;...
    perCaUMin, perCaUMax;...
    perCaOHMin, perCaOHMax;...
    perCaCOMin, perCaCOMax];

%% perovskite, MgFe

% 0.01 to 0.5 for Sm
% 0.01 to 0.02 for Nd
% 0.01 to 1 for U
% 0.01 to 0.1 for Th

perMgFeSmMin = 0.01;
perMgFeSmMax = 0.5;

perMgFeNdMin = 0.01;
perMgFeNdMax = 0.1;

perMgFeUMin = 0.01;
perMgFeUMax = 1.1;

perMgFeThMin = 0.01;
perMgFeThMax = 0.2;

perMgFeOHMin = 0.00001;
perMgFeOHMax = 0.001;

perMgFeCOMin = 0.00005;
perMgFeCOMax = 0.005;

Kd_p_MgFe_minmax = [perMgFeSmMin, perMgFeSmMax;...
    perMgFeNdMin, perMgFeNdMax;...
    perMgFeThMin, perMgFeThMax ;...
    perMgFeUMin, perMgFeUMax;...
    perMgFeOHMin, perMgFeOHMax;...
    perMgFeCOMin, perMgFeCOMax];

%% post-perovskite

%% loop

% final matrix:
% run number - 1
% del-LJ values - N - 1, N - 2

%printableRunName = cell(numRuns, 1); %zeros(numRuns, 1);
runInfo = zeros(numRuns, (kdNum + 3 + 4)); % + 4 for U & Th

for loop = 1:1:numRuns
    %runName = cellstr(strcat('run', num2str(i)));
    %printableRunName(i, 1) = runName;
    runInfo(loop, 1) = loop;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Magnesiowustite
    
%     % Sm
%     %outValues(i, 1)
%     runInfo(loop, 2) = Kd_m_minmax(1,1) +...
%         (Kd_m_minmax(1,2) - Kd_m_minmax(1,1)).*rand(1);
%     KSm_m = runInfo(loop, 2);
%     
%     % Nd
%     %outValues(i, 2)
%     runInfo(loop, 3) = Kd_m_minmax(2,1) +...
%         (Kd_m_minmax(2,2) - Kd_m_minmax(2,1)).*rand(1);
%     KNd_m = runInfo(loop, 3);
%     
%     % Th
%     %outValues(i, 3)
%     runInfo(loop, 4) = Kd_m_minmax(3,1) +...
%         (Kd_m_minmax(3,2) - Kd_m_minmax(3,1)).*rand(1);
%     KTh_m = runInfo(loop, 4);
%     
%     % U
%     %outValues(i, 4)
%     runInfo(loop, 5) = Kd_m_minmax(4,1) +...
%         (Kd_m_minmax(4,2) - Kd_m_minmax(4,1)).*rand(1);
%     KU_m = runInfo(loop, 5);
%     
%     % OH
%     %outValues(i, 5)
%     runInfo(loop, 6) = Kd_m_minmax(5,1) +...
%         (Kd_m_minmax(5,2) - Kd_m_minmax(5,1)).*rand(1);
%     KOH_m = runInfo(loop, 6);
%     
%     % CO
%     %outValues(i, 6)
%     runInfo(loop, 7) = Kd_m_minmax(6,1) +...
%         (Kd_m_minmax(6,2) - Kd_m_minmax(6,1)).*rand(1);
%     KCO_m = runInfo(loop, 7);
%     
%     Kd_m = [KSm_m; KNd_m; KTh_m; KU_m; KOH_m; KCO_m];
    Kd_m = [];
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Ca-Perovskite
    runInfo(loop, 8) = Kd_p_Ca_minmax(1,1) +...
        (Kd_p_Ca_minmax(1,2) - Kd_p_Ca_minmax(1,1)).*rand(1);
    KSm_p_Ca = runInfo(loop, 8);
    
    % Nd
    %outValues(loop, 2)
    runInfo(loop, 9) = Kd_p_Ca_minmax(2,1) +...
        (Kd_p_Ca_minmax(2,2) - Kd_p_Ca_minmax(2,1)).*rand(1);
    KNd_p_Ca = runInfo(loop, 9);
    
    % Th
    %outValues(loop, 3)
    runInfo(loop, 10) = Kd_p_Ca_minmax(3,1) +...
        (Kd_p_Ca_minmax(3,2) - Kd_p_Ca_minmax(3,1)).*rand(1);
    KTh_p_Ca = runInfo(loop, 10);
    
    % U
    %outValues(loop, 4)
    runInfo(loop, 11) = Kd_p_Ca_minmax(4,1) +...
        (Kd_p_Ca_minmax(4,2) - Kd_p_Ca_minmax(4,1)).*rand(1);
    KU_p_Ca = runInfo(loop, 12);
    
    % OH
    %outValues(loop, 5)
    runInfo(loop, 12) = Kd_p_Ca_minmax(5,1) +...
        (Kd_p_Ca_minmax(5,2) - Kd_p_Ca_minmax(5,1)).*rand(1);
    KOH_p_Ca = runInfo(loop, 12);
    
    % CO
    %outValues(loop, 6)
    runInfo(loop, 13) = Kd_p_Ca_minmax(6,1) +...
        (Kd_p_Ca_minmax(6,2) - Kd_p_Ca_minmax(6,1)).*rand(1);
    KCO_p_Ca = runInfo(loop, 13);
    
    Kd_p_Ca = [KSm_p_Ca; KNd_p_Ca; KTh_p_Ca; KU_p_Ca; KOH_p_Ca; KCO_p_Ca];  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MgFe-Perovskite
    runInfo(loop, 14) = Kd_p_MgFe_minmax(1,1) +...
        (Kd_p_MgFe_minmax(1,2) - Kd_p_MgFe_minmax(1,1)).*rand(1);
    KSm_p_MgFe = runInfo(loop, 14);
    
    % Nd
    %outValues(loop, 2)
    runInfo(loop, 15) = Kd_p_MgFe_minmax(2,1) +...
        (Kd_p_MgFe_minmax(2,2) - Kd_p_MgFe_minmax(2,1)).*rand(1);
    KNd_p_MgFe = runInfo(loop, 15);
    
    % Th
    %outValues(loop, 3)
    runInfo(loop, 16) = Kd_p_MgFe_minmax(3,1) +...
        (Kd_p_MgFe_minmax(3,2) - Kd_p_MgFe_minmax(3,1)).*rand(1);
    KTh_p_MgFe = runInfo(loop, 16);
    
    % U
    %outValues(loop, 4)
    runInfo(loop, 17) = Kd_p_MgFe_minmax(4,1) +...
        (Kd_p_MgFe_minmax(4,2) - Kd_p_MgFe_minmax(4,1)).*rand(1);
    KU_p_MgFe = runInfo(loop, 17);
    
    % OH
    %outValues(loop, 5)
    runInfo(loop, 18) = Kd_p_MgFe_minmax(5,1) +...
        (Kd_p_MgFe_minmax(5,2) - Kd_p_MgFe_minmax(5,1)).*rand(1);
    KOH_p_MgFe = runInfo(loop, 18);
    
    % CO
    %outValues(loop, 6)
    runInfo(loop, 19) = Kd_p_MgFe_minmax(6,1) +...
        (Kd_p_MgFe_minmax(6,2) - Kd_p_MgFe_minmax(6,1)).*rand(1);
    KCO_p_MgFe = runInfo(loop, 19);
    
    Kd_p_MgFe = [KSm_p_MgFe; KNd_p_MgFe; KTh_p_MgFe; KU_p_MgFe; ...
        KOH_p_MgFe; KCO_p_MgFe];  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %     delEDR = 0;
    %     delEER = -60;
    
    
    all_oceans

    %% EDR & EER
    runInfo(loop, 20) = EDR_del_LJ_print;
    runInfo(loop, 21) = EER_del_LJ_print;

    %% U & Th
    runInfo(loop, 22) = EERUfracoftotal;
    runInfo(loop, 23) = EERThfracoftotal;
    runInfo(loop, 24) = EERUfracoftotalwliq;
    runInfo(loop, 25) = EERThfracoftotalwliq;

end

toc

dlmwrite('runInfo.dat', runInfo, 'precision', '%.6f', ...
         'newline', 'pc')
