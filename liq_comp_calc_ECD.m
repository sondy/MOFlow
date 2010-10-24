%%
% liq_comp_calc_ECD.m
% October 23, 2010
% Alessondra Springmann
% inputs:
% liq_comp, mass_liquid, mass_solidified, mass_this_shell, new_mass_liquid
% solid, liquid

function new_liq = liq_comp_calc_ECD(...
    j, ...
    delr, ...
    liquid, ... 
    Mantlevolume, ...
    mass_liquid, ...
    mass_solidified, ...
    mass_this_shell, ...
    new_mass_liquid, ...
    solid)

liq_term = liquid(j-1,:) .* (mass_liquid./new_mass_liquid);

sol_term = (solid(j,:) .*...
    ((mass_solidified - mass_this_shell)./new_mass_liquid));

new_liq = liq_term.*(1 + (Mantlevolume + delr(j)./Mantlevolume)) ...
    + liq_term ...
    - sol_term;