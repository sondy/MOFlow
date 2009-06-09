function VP = BirchMurnsolid(Kot, Kpt, Vo, P)
% Birch-Murnaghan solution for solid density with depth
% calculates VP at given P based on Vo (which is V(t) at one bar)
% technique from Fei et al (1990) etc (but note that many paper's equations 
% are wrong)
% P is pressure in GPa but Birch-Murn calcs in Pa

Pfunc = @(Vratio)(3/2)*Kot*((Vratio)^(7/3) - ...
    (Vratio)^(5/3))*(1 - (3/4)*(4 - Kpt)*((Vratio)^(2/3) - 1)) - P;

range = [1 4];
% range = 1.3;

[Vpressure, FVAL, EXITFLAG] = fzero(Pfunc, range);

global xrecord xrecordn
xrecordn = xrecordn + 1;
xrecord(xrecordn) = Vpressure;

if (EXITFLAG < 0)
    fprintf('\n');
    fprintf('fzero failed --- Vratio guess is less than zero \n');
    fprintf('\n');
    PrintCaller
    fprintf('\n');
    fprintf('\n');
end

VP = Vo/Vpressure;

end

%VP = Vo/(fzero(Pfunc,1.2));

% Pfunc is just the birch-murn equation with Vratio in the place of Vo/VP
% and the entire equation rearranged such that Pfunc is zero at the correct
% volume ratio.
% If Vratio = Vo/VP, then VP = Vo/Vratio
% [0.95 4] for the initial guess seems to be a good one for the entire 
% volume of the Earth's mantle for perovskite
% fzero simply finds the zero point of the Pfunc function at around 1.2 and
% returns Vratio

% If speed ever becomes a problem, calculate the inverse of this equation
% symbolically then use it here

% This doesn't seem to work for other minerals in the rest of the mantle, as
% it begins returning imaginary results!  Hrm.  I wonder how to best generate
% guesses for other minerals; I guess looking at how Vo/VP behaves for other
% minerals would be a good idea.

% CalcPressure = (3/2)*Kot*((Vo/VP)^(7/3) - (Vo/VP)^(5/3))*(1 - (3/4)*(4 -
% Kpt)*((Vo/VP)^(2/3) - 1));
