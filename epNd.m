% epNd.m
% 1/4/2011
% Alessondra Springmann
% The goal: calculate epsilon 142/144 Nd for the model
% This script can run independently of MOFlowCMB --- not any more!

%chonEp142Nd = 1.141901; % from Amelin, 2004
%chonEp142Nd = 1.149498;

% Carlson spreadsheet
Nd142Nd144 = zeros(1, 4567);
Nd142Nd144(1) = 1.1414980; %B2
B2 = Nd142Nd144(1);

Sm144Nd144 = 0.0402; %E2
E2 = Sm144Nd144;

Sm147Nd144 = 0.19600; %F2
F2 = Sm147Nd144;

standard = 1.14184; %G2
G2 = standard;

mysteryRatio = 0.51263; %G3
G3 = mysteryRatio;

delChon = -18.0; %H2
H2 = delChon;

chon142144 = 1.141819447; %I2
I2 = chon142144;

Nd142144Nd = G3-F2*(exp(0.00000000000654*4567000000)-1); %K2
K2 = Nd142144Nd;

mysteryRatio2 = 0.51263; %K3
K3 = mysteryRatio2;

% from data tables on Nd
abunNd142 = 0.272; % number ratio abundance "atom percent abundance"

%abunNd144 = 0.238; % number ratio abundance

atomicRatio147 = 0.20503;

abunNd142fromSm = 0.008; % number ratio abundance from Scott, 2007
C2 = abunNd142fromSm;
% this is from Scott (2007): ratio of 146/144 Sm = 8 * 10^-3
% [146] + [144] = 1; [146]/[144] = 8e-3
% so [146] has to be 0.008; [144] = 0.992
% all of 146 Sm has decayed to 142 Nd

decayConstant = 6.73e-09;
D2 = decayConstant;

massNd142 = 141.907719; % atomic mass 142-neodymium

massNd144 = 143.910083; % atomic mass 144-neodymium

massSm = 150.36;

massNd = 144.24;

abunSm147 = 0.1499;

abunNd144 = 0.238;

%% EDR CONSTANTS

EDR_wtpercentSm = 0.0000049307;
EDR_molesSm = EDR_wtpercentSm/massSm;
EDR_moles147Sm = EDR_molesSm*abunSm147;

EDR_wtpercentNd = 0.000008722;
EDR_molesNd = EDR_wtpercentNd/massNd;
EDR_moles144Nd = EDR_molesNd*abunNd144;

EDR_147Sm144Nd = EDR_moles147Sm/EDR_moles144Nd;
F3 = EDR_147Sm144Nd;
E3 = F3*abunNd142fromSm;

%% EER CONSTANTS

EER_wtpercentSm = 0.00013856;
EER_molesSm = EER_wtpercentSm/massSm;
EER_moles147Sm = EER_molesSm*abunSm147;

EER_wtpercentNd = 0.00050067;
EER_molesNd = EER_wtpercentNd/massNd;
EER_moles144Nd = EER_molesNd*abunNd144;

EER_147Sm144Nd = EER_moles147Sm/EER_moles144Nd;
F4 = EER_147Sm144Nd;
E4 = F4*abunNd142fromSm;

%% time

% COLUMN A
T = 1:1:4567;

maxsize = length(T);

% Rest of row 7

% COLUMN C
del_LJ = zeros(1, 4567);
del_LJ(1) = ((B2/G2) - 1)*1000000;

% COLUMN D
Sm146Sm144 = zeros(1, 4567);
Sm146Sm144(1) = C2*exp(-D2*T(1)*1000000);

% COLUMN E
Nd143Nd144 = zeros(1, 4567);
Nd143Nd144(1) = K2;

% COLUMN F
Sm147Nd144 = zeros(1, 4567);
Sm147Nd144(1) = F2/(exp(-1*0.00000000000654*(4567000000-T(1)*1000000)));

%% EDR
% COLUMN G
EDR_Nd142Nd144 = zeros(1, 4567);
EDR_Nd142Nd144(1) = B2;

% COLUMN H
EDR_del_LJ = zeros(1, 4567);
EDR_del_LJ(1) = (EDR_Nd142Nd144(1)/G2-1)*1000000; % G7 = EDR_Nd142Nd144;

% COLUMN I
% EDR_143144Nd = zeros(1, 4567);
% EDR_143144Nd(1) = K2;

%% EER

% COLUMN L
EER_Nd142Nd144 = zeros(1, 4567);
EER_Nd142Nd144(1) = B2;

% COLUMN M
EER_del_LJ = zeros(1, 4567);
EER_del_LJ(1) = (EER_Nd142Nd144(1)/G2-1)*1000000; % L7 = EDR_Nd142Nd144;

%% THE BIG LOOP, with liquids
for j = 2:1:maxsize
    % COLUMN B
    Nd142Nd144(j) = B2+((E2*C2)*((1-exp(-D2*T(j)*1000000))));
    
    % COLUMN C
    del_LJ(j) = (Nd142Nd144(j)/G2-1)*1000000;
    
    % COLUMN D
    Sm146Sm144(j) = C2*exp(-D2*T(j)*1000000);
    
    % COLUMN E
    Nd143Nd144(j) = G3-F2*(exp(0.00000000000654*(4567-T(j))*1000000)-1);
    
    % COLUMN F
    Sm147Nd144(j) = F2/(exp(-1*0.00000000000654*(4567000000-T(j)*1000000)));
    
    %% EDR
    % COLUMN G
    EDR_Nd142Nd144(j) = EDR_Nd142Nd144(j-1) + ((E3*Sm146Sm144(j-1))*((1-exp(-D2*1000000))));
    % =G7+((E$3*D7)*((1-EXP(-D$2*(A8-A7)*1000000))))
    
    % COLUMN H
    EDR_del_LJ(j) = (EDR_Nd142Nd144(j)/G2-1)*1000000;
    
    %% EER
    
    % COLUMN L
    EER_Nd142Nd144(j) = EER_Nd142Nd144(j-1) + ((E4*Sm146Sm144(j-1))*((1-exp(-D2*1000000))));
    
    % COLUMN M
    EER_del_LJ(j) = (EER_Nd142Nd144(j)/G2-1)*1000000;
end


%% THE BIG PRINT

fprintf('The final time is %3g Myr. CHECK \n', T(maxsize)) % CHECK
fprintf('The final 142Nd/144Nd ratio is %2.8g. CHECK \n', Nd142Nd144(maxsize))
fprintf('The final del-LJ is %2.3g. CHECK \n', del_LJ(maxsize))
fprintf('The final 146Sm/144Sm ratio is %2.6g. CHECK \n', 0); %Sm146Sm144(maxsize)) % effectively zero.
fprintf('The final 143Nd/144Nd ratio is %2.6g. CHECK \n', Nd143Nd144(maxsize))
fprintf('The final 147Sm/144Nd ratio is %2.6g. CHECK \n', Sm147Nd144(maxsize))
fprintf('The final EDR 142Nd/144Nd ratio is %2.8g. \n', EDR_Nd142Nd144(maxsize))
fprintf('The final EDR del-LJ ratio is %2.6g. \n', EDR_del_LJ(maxsize))
fprintf('The final EER 142Nd/144Nd ratio is %2.8g. \n', EER_Nd142Nd144(maxsize))
fprintf('The final EER del-LJ ratio is %2.6g. \n', EER_del_LJ(maxsize))
fprintf('\n')