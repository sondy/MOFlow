% initialatmospherewater.m

disp('In initialatmospherewater')

% Given the initial bulk composition of the magma ocean in terms of water, 
% iterates to partition that quantity into the magma ocean and the
% atmosphere for an equilibrium saturation to begin the model from

%%%calculate initial conditions
% H2Oliquid(1) is variable name from FLOW program
% H2Oliq is variable that changes within this program
% same for HPatm vs HPa and HMatm vs HMa

Mwater = (H2Oliquid(1)/100)*Mantlemass;   % total mass of water and CO2 in kg in the system

H2Oliqinit = 0.1*H2Oliquid(1);  % first guess at wt% in MO (as opposed to atm) MUST be less then actual

HPa = [((H2Oliqinit - 0.30)/(2.08e-6*100))^(1/0.52)];   % leads to first guess for pressure in atm
        if isreal(HPa) == 0; HPa = 0; end                   % in case initial value is below 0.3
        if H2Oliqinit < 0.3; HPa = 0; end

HM = (HPa*4*pi*R^2)/g + (H2Oliqinit/100)*Mantlemass;    % first guess of mass in atmosphere + MO

Hadjust = H2Oliqinit*0.03; % amount to add or subtract from H2Oliq if error is too large, in wt%

Herror = Mwater*0.01; % maximum allowable error in kg, ends iterative solution loop

%% First for water
while abs(HM - Mwater) >= Herror % if error is too large, iterate
        if (HM - Mwater) > 0 % if HM is too high, H2Oliqinit is too high
          H2Oliqinit = H2Oliqinit - Hadjust;
        elseif (HM - Mwater) < 0
          H2Oliqinit = H2Oliqinit + Hadjust;
        end
    HPa = [((H2Oliqinit - 0.30)/(2.08e-6*100))^(1/0.52)];   % water partial pressure in Pa
        if isreal(HPa) == 0; HPa = 0; end                   % in case initial value is below 0.3
        if H2Oliqinit < 0.3; HPa = 0; end
    HM = (HPa*4*pi*R^2)/g + (H2Oliqinit/100)*Mantlemass;          % next guess of mass in atmosphere + MO
end

H2Oliquid2(1) = H2Oliqinit;
HPatm(1) = HPa;
