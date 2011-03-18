% epNd.m
% 1/4/2011
% Alessondra Springmann
% The goal: calculate epsilon 142/144 Nd for the model
% This script can run independently of MOFlowCMB --- not any more!

%chonEp142Nd = 1.141901; % from Amelin, 2004
%chonEp142Nd = 1.149498;

% Carlson spreadsheet

%clear all;

%% Time

% COLUMN A
T = 0:1:4567;

maxsize = length(T);

Myr = 1000000;

%% Important Constants

atomicRatio147 = 0.20503;

% this is from Scott (2007): ratio of 146/144 Sm = 8 * 10^-3
% [146] + [144] = 1; [146]/[144] = 8e-3
% so [146] has to be 0.008; [144] = 0.992
% all of 146 Sm has decayed to 142 Nd
abunNd142fromSm = 0.008; % number ratio abundance from Scott, 2007
C2 = abunNd142fromSm;

Sm147Nd144 = 0.19600; %F2
F2 = Sm147Nd144;

Sm144Nd144 = F2*atomicRatio147; %E2
E2 = Sm144Nd144;

standard = 1.14184; %G2
G2 = standard;

mysteryRatio = 0.51263; %G3
G3 = mysteryRatio;

delChon = -18.0; %H2
H2 = delChon;

chon142144 = (H2./Myr+1)*G2; %I2
I2 = chon142144;

mysteryConstant = 0.00000000000654;

Nd143Nd144 = G3-F2*(exp(mysteryConstant*maxsize*Myr)-1); %K2
K2 = Nd143Nd144;

mysteryRatio2 = 0.51263; %K3
K3 = mysteryRatio2;

Nd142Nd144 = zeros(1, 4567);
Nd142Nd144(1) = I2-E2.*C2; %B2
%B2 = Nd142Nd144(1);
B2 = Nd142Nd144(1);

% from data tables on Nd
abunNd142 = 0.272; % number ratio abundance "atom percent abundance"

%abunNd144 = 0.238; % number ratio abundance

decayConstant = 6.73e-09;
D2 = decayConstant;

massNd142 = 141.907719; % atomic mass 142-neodymium

massNd144 = 143.910083; % atomic mass 144-neodymium

massSm = 150.36;

massNd = 144.24;

abunSm147 = 0.1499;

abunNd144 = 0.238;

%% EDR CONSTANTS

EDR_wtpercentSm = avgSmEDR; %0.0000049307;
EDR_molesSm = EDR_wtpercentSm/massSm;
EDR_moles147Sm = EDR_molesSm*abunSm147;

EDR_wtpercentNd = avgNdEDR; %0.000008722;
EDR_molesNd = EDR_wtpercentNd/massNd;
EDR_moles144Nd = EDR_molesNd*abunNd144;

EDR_147Sm144Nd = EDR_moles147Sm/EDR_moles144Nd;
F3 = EDR_147Sm144Nd;
E3 = F3*atomicRatio147;

%% EER CONSTANTS

EER_wtpercentSm = avgSmEERwLiq; %0.00013856;
EER_molesSm = EER_wtpercentSm/massSm;
EER_moles147Sm = EER_molesSm*abunSm147;

EER_wtpercentNd = avgNdEERwLiq; %0.00050067;
EER_molesNd = EER_wtpercentNd/massNd;
EER_moles144Nd = EER_molesNd*abunNd144;

EER_147Sm144Nd = EER_moles147Sm/EER_moles144Nd;
F4 = EER_147Sm144Nd;
E4 = F4*atomicRatio147;


%% Initial values; rest of row 7

% COLUMN C
del_LJ = zeros(1, maxsize);
del_LJ(1) = ((B2/G2) - 1)*Myr;

% COLUMN D
Sm146Sm144 = zeros(1, maxsize);
Sm146Sm144(1) = C2*exp(-D2*T(1)*Myr);

% COLUMN E
Nd143Nd144 = zeros(1, maxsize);
Nd143Nd144(1) = K2;

% COLUMN F
Sm147Nd144 = zeros(1, maxsize);
Sm147Nd144(1) = F2/(exp(-1*mysteryConstant*(4567000000-T(1)*Myr)));

%% EDR
% COLUMN G
EDR_Nd142Nd144 = zeros(1, maxsize);
EDR_Nd142Nd144(1) = B2;

% COLUMN H
EDR_del_LJ = zeros(1, maxsize);
EDR_del_LJ(1) = (EDR_Nd142Nd144(1)/G2-1)*Myr; % G7 = EDR_Nd142Nd144;

% COLUMN I
% EDR_143144Nd = zeros(1, 4567);
% EDR_143144Nd(1) = K2;

%% EER

% COLUMN L
EER_Nd142Nd144 = zeros(1, maxsize);
EER_Nd142Nd144(1) = B2;

% COLUMN M
EER_del_LJ = zeros(1, maxsize);
EER_del_LJ(1) = (EER_Nd142Nd144(1)/G2-1)*Myr; % L7 = EDR_Nd142Nd144;

%% initial values
% fprintf('\n')
% fprintf('\n')
% fprintf('The initial time is %3g Myr.  \n', T(1)) % CHECK
% fprintf('The initial 142Nd/144Nd ratio is %2.9g.  \n', Nd142Nd144(1))
% fprintf('The initial del-LJ is %2.4g.  \n', del_LJ(1))
% fprintf('The initial 146Sm/144Sm ratio is %2.9g.  \n', Sm146Sm144(1)); %Sm146Sm144(maxsize)) % effectively zero.
% fprintf('The initial 143Nd/144Nd ratio is %2.9g.  \n', Nd143Nd144(1))
% fprintf('The initial 147Sm/144Nd ratio is %2.9g.  \n', Sm147Nd144(1))
% fprintf('The initial EDR 142Nd/144Nd ratio is %2.9g. \n', EDR_Nd142Nd144(1))
% fprintf('The initial EDR del-LJ ratio is %2.4g. \n', EDR_del_LJ(1))
% fprintf('The initial EER 142Nd/144Nd ratio is %2.9g. \n', EER_Nd142Nd144(1))
% fprintf('The initial EER del-LJ ratio is %2.4g. \n', EER_del_LJ(1))
% fprintf('\n')
% 
% check = [];

%% THE BIG LOOP, with liquids
for j = 2:1:maxsize
    % COLUMN B
    Nd142Nd144(j) = B2+((E2*C2)*((1-exp(-D2*T(j)*Myr))));
    
    % COLUMN C
    del_LJ(j) = (Nd142Nd144(j)/G2-1)*Myr;
    
    % COLUMN D
    Sm146Sm144(j) = C2*exp(-D2*T(j)*Myr);
    
    % COLUMN E
    Nd143Nd144(j) = G3-F2*(exp(mysteryConstant*(4567-T(j))*Myr)-1);
    
    % COLUMN F
    Sm147Nd144(j) = F2/exp(-1*mysteryConstant*(4567-T(j))*Myr);
    
    %% EDR
    % COLUMN G
    EDR_Nd142Nd144(j) = EDR_Nd142Nd144(j-1) + E3.*Sm146Sm144(j-1).*(1-exp(-D2.*(T(j)-T(j-1)).*Myr));

    % COLUMN H
    EDR_del_LJ(j) = ((EDR_Nd142Nd144(j)/G2)-1)*Myr;
    
    %% EER
    
    % COLUMN L
    EER_Nd142Nd144(j) = EER_Nd142Nd144(j-1) + ((E4.*Sm146Sm144(j-1)).*((1-exp(-D2.*(T(j)-T(j-1)).*Myr))));
    
    % COLUMN M
    EER_del_LJ(j) = ((EER_Nd142Nd144(j)/G2)-1)*Myr;
end


%% THE BIG PRINT

% fprintf('The final time is %3g Myr.\n', T(maxsize)) % CHECK
% fprintf('The final 142Nd/144Nd ratio is %2.8g.\n', Nd142Nd144(maxsize))
% fprintf('The final del-LJ is %2.4g.\n', del_LJ(maxsize))
% fprintf('The final 146Sm/144Sm ratio is %0.1g.\n', Sm146Sm144(maxsize)); %Sm146Sm144(maxsize)) % effectively zero.
% fprintf('The final 143Nd/144Nd ratio is %2.6g.\n', Nd143Nd144(maxsize))
% fprintf('The final 147Sm/144Nd ratio is %2.6g.\n', Sm147Nd144(maxsize))
% fprintf('The final EDR 142Nd/144Nd ratio is %2.8g. \n', EDR_Nd142Nd144(maxsize))
fprintf('\n \n')
fprintf('The final EDR del-LJ ratio is %2.4g. \n', EDR_del_LJ(maxsize))
% fprintf('The final EER 142Nd/144Nd ratio is %2.8g. \n', EER_Nd142Nd144(maxsize))
% fprintf('The final EER del-LJ ratio is %2.4g. \n', EER_del_LJ(maxsize))
fprintf('\n')
