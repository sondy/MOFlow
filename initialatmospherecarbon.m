% initialatmospherecarbon.m

disp('In initialatmospherecarbon')

% Given the initial bulk composition of the magma ocean in terms of CO2, 
% iterates to partition that quantity into the magma ocean and the
% atmosphere for an equilibrium saturation to begin the model from

%%%calculate initial conditions
% H2Oliquid(1) is variable name from FLOW program
% H2Oliq is variable that changes within this program
% same for HPatm vs HPa and HMatm vs HMa

Mcarbon = (CO2liquid(1)/100)*Mantlemass;  % this is what needs to be matched by calculations here

CO2liqinit = 0.1*CO2liquid(1);

CPa = [((CO2liqinit - 0.05)/(2.08e-6*100))^(1/0.45)];   % converts mass percent to pressure in Pa
        if isreal(CPa) == 0; CPa = 0; end                   % in case initial value is below 0.05
        if CO2liqinit < 0.05; CPa = 0; end

CM = (CPa*4*pi*R^2)/g + (CO2liqinit/100)*Mantlemass;   % translates P in Pa to Matm in kg

Cadjust = CO2liqinit*0.03;

Cerror = Mcarbon*0.01;

% for CO2
while abs(CM - Mcarbon) >= Cerror % if error is too large, iterate
        if (CM - Mcarbon) > 0 % if CM is too high, CO2liqinit is too high
          CO2liqinit = CO2liqinit - Cadjust;
        elseif (CM - Mcarbon) < 0
          CO2liqinit = CO2liqinit + Cadjust;
        end
    CPa = [((CO2liqinit - 0.05)/(2.08e-6*100))^(1/0.45)];   % CO2 partial pressure in Pa
        if isreal(CPa) == 0; CPa = 0; end                   % in case initial value is below 0.05
        if CO2liqinit < 0.05; CPa = 0; end
    CM = (CPa*4*pi*R^2)/g + (CO2liqinit/100)*Mantlemass;          % leads to first guess of mass in atmosphere + MO
end

CO2liquid2(1) = CO2liqinit;
CPatm(1) = CPa;