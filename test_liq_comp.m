liq_comp = [1 2 3 4 5 6 7 8 9 10 11];

maxstep = 997;

liquid_composition = zeros(maxstep, length(liq_comp));

for j = 1:1:maxstep;
   liquid_composition(j, :) = liq_comp;
end

display(liquid_composition)