% monteLoop.m
% Alessondra Springmann, 3/19/2011
% Monte Carlo simulations for MOFlow, specifically, changing the partition
% coefficients for perovskite, post-perovskite, and magnesiowustite

% decide how many runs will occur
% create a range of variables for each run
% write run number, input variables into a file
% calculate output values, write into a file

numRuns = 10;

numPar = 6;

numMin = 3;

kdNum = numPar*numMin; % number of partition coefficients * number of minerals

kdMatrix = cell(numRuns, kdNum);

%% magnesiowustite
magSmNdUMin = 0.0001;
magSmNdUMax = 0.1;

magThMin = 0.001;
magThMax = 1.0;

magOHMin = 0.0001;
magOHMax = 0.1;

magCOMin = 0.00001;
magCOMax = 0.01;

Kd_m_minmax = [magSmNdUMin, magSmNdUMax;...
    magSmNdUMin, magSmNdUMin;...
    magThMin, magThMax ;...
    magSmNdUMin, magSmNdUMax;...
    magOHMin, magOHMax;...
    magCOMin, magCOMax];

%% perovskite, Ca

% 1 - 20 for Sm
% 5 - 20 for Nd
% 0.01 - 20 for U
% 1 - 20 for Th

perCaSmMin = 0.9;
perCaSmMax = 20;

perCaNdMin = 4;
perCaNdMax = 20;

perCaUMin = 0.009;
perCaUMax = 20;

perCaThMin = 0.9;
perCaThMax = 20;

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

perMgFeSmMin = 0.009;
perMgFeSmMax = 0.6;

perMgFeNdMin = 0.009;
perMgFeNdMax = 0.03;

perMgFeUMin = 0.009;
perMgFeUMax = 1.1;

perMgFeThMin = 0.009;
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
runInfo = zeros(numRuns, (kdNum + 3));

for i = 1:1:numRuns
    %runName = cellstr(strcat('run', num2str(i)));
    %printableRunName(i, 1) = runName;
    runInfo(i, 1) = i;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Magnesiowustite
    
    % Sm
    %outValues(i, 1)
    runInfo(i, 2) = Kd_m_minmax(1,1) +...
        (Kd_m_minmax(1,2) - Kd_m_minmax(1,1)).*rand(1);
    KSm_m = runInfo(i, 2);
    
    % Nd
    %outValues(i, 2)
    runInfo(i, 3) = Kd_m_minmax(2,1) +...
        (Kd_m_minmax(2,2) - Kd_m_minmax(2,1)).*rand(1);
    KNd_m = runInfo(i, 3);
    
    % Th
    %outValues(i, 3)
    runInfo(i, 4) = Kd_m_minmax(3,1) +...
        (Kd_m_minmax(3,2) - Kd_m_minmax(3,1)).*rand(1);
    KTh_m = runInfo(i, 4);
    
    % U
    %outValues(i, 4)
    runInfo(i, 5) = Kd_m_minmax(4,1) +...
        (Kd_m_minmax(4,2) - Kd_m_minmax(4,1)).*rand(1);
    KU_m = runInfo(i, 5);
    
    % OH
    %outValues(i, 5)
    runInfo(i, 6) = Kd_m_minmax(5,1) +...
        (Kd_m_minmax(5,2) - Kd_m_minmax(5,1)).*rand(1);
    KOH_m = runInfo(i, 6);
    
    % CO
    %outValues(i, 6)
    runInfo(i, 7) = Kd_m_minmax(6,1) +...
        (Kd_m_minmax(6,2) - Kd_m_minmax(6,1)).*rand(1);
    KCO_m = runInfo(i, 7);
    
    Kd_m = [KSm_m; KNd_m; KTh_m; KU_m; KOH_m; KCO_m];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Ca-Perovskite
    runInfo(i, 8) = Kd_p_Ca_minmax(1,1) +...
        (Kd_p_Ca_minmax(1,2) - Kd_p_Ca_minmax(1,1)).*rand(1);
    KSm_p_Ca = runInfo(i, 2);
    
    % Nd
    %outValues(i, 2)
    runInfo(i, 9) = Kd_p_Ca_minmax(2,1) +...
        (Kd_p_Ca_minmax(2,2) - Kd_p_Ca_minmax(2,1)).*rand(1);
    KNd_p_Ca = runInfo(i, 3);
    
    % Th
    %outValues(i, 3)
    runInfo(i, 10) = Kd_p_Ca_minmax(3,1) +...
        (Kd_p_Ca_minmax(3,2) - Kd_p_Ca_minmax(3,1)).*rand(1);
    KTh_p_Ca = runInfo(i, 4);
    
    % U
    %outValues(i, 4)
    runInfo(i, 11) = Kd_p_Ca_minmax(4,1) +...
        (Kd_p_Ca_minmax(4,2) - Kd_p_Ca_minmax(4,1)).*rand(1);
    KU_p_Ca = runInfo(i, 5);
    
    % OH
    %outValues(i, 5)
    runInfo(i, 12) = Kd_p_Ca_minmax(5,1) +...
        (Kd_p_Ca_minmax(5,2) - Kd_p_Ca_minmax(5,1)).*rand(1);
    KOH_p_Ca = runInfo(i, 6);
    
    % CO
    %outValues(i, 6)
    runInfo(i, 13) = Kd_p_Ca_minmax(6,1) +...
        (Kd_p_Ca_minmax(6,2) - Kd_p_Ca_minmax(6,1)).*rand(1);
    KCO_p_Ca = runInfo(i, 7);
    
    Kd_p_Ca = [KSm_p_Ca; KNd_p_Ca; KTh_p_Ca; KU_p_Ca; KOH_p_Ca; KCO_p_Ca];  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MgFe-Perovskite
    runInfo(i, 14) = Kd_p_MgFe_minmax(1,1) +...
        (Kd_p_MgFe_minmax(1,2) - Kd_p_MgFe_minmax(1,1)).*rand(1);
    KSm_p_MgFe = runInfo(i, 2);
    
    % Nd
    %outValues(i, 2)
    runInfo(i, 15) = Kd_p_MgFe_minmax(2,1) +...
        (Kd_p_MgFe_minmax(2,2) - Kd_p_MgFe_minmax(2,1)).*rand(1);
    KNd_p_MgFe = runInfo(i, 3);
    
    % Th
    %outValues(i, 3)
    runInfo(i, 16) = Kd_p_MgFe_minmax(3,1) +...
        (Kd_p_MgFe_minmax(3,2) - Kd_p_MgFe_minmax(3,1)).*rand(1);
    KTh_p_MgFe = runInfo(i, 4);
    
    % U
    %outValues(i, 4)
    runInfo(i, 17) = Kd_p_MgFe_minmax(4,1) +...
        (Kd_p_MgFe_minmax(4,2) - Kd_p_MgFe_minmax(4,1)).*rand(1);
    KU_p_MgFe = runInfo(i, 5);
    
    % OH
    %outValues(i, 5)
    runInfo(i, 18) = Kd_p_MgFe_minmax(5,1) +...
        (Kd_p_MgFe_minmax(5,2) - Kd_p_MgFe_minmax(5,1)).*rand(1);
    KOH_p_MgFe = runInfo(i, 6);
    
    % CO
    %outValues(i, 6)
    runInfo(i, 19) = Kd_p_MgFe_minmax(6,1) +...
        (Kd_p_MgFe_minmax(6,2) - Kd_p_MgFe_minmax(6,1)).*rand(1);
    KCO_p_MgFe = runInfo(i, 7);
    
    Kd_p_MgFe = [KSm_p_MgFe; KNd_p_MgFe; KTh_p_MgFe; KU_p_MgFe; ...
        KOH_p_MgFe; KCO_p_MgFe];  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %     delEDR = 0;
    %     delEER = -60;
    
    all_oceans
    
    runInfo(i, 20) = EDR_del_LJ;
    runInfo(i, 21) = EER_del_LJ;
    
end