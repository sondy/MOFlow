%%
% liq_comp_calc_HMH.m
% October 23, 2010
% Alessondra Springmann
% inputs:
% liq_comp, mass_liquid, mass_solidified, mass_this_shell, new_mass_liquid
% solid, liquid

function new_liq = liq_comp_calc_HMH(j, liquid, mass_liquid,...
    mass_this_shell, new_mass_liquid, solid)

new_liq = (liquid(j-1,:) * mass_liquid...
    - mass_this_shell * solid(j,:))/new_mass_liquid;