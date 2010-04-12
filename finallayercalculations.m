% finallayercalculations.m

% common final calculations after each increment of solid fractionation

%disp('Final Layer Calculations')

mass_this_shell = density*(Mantlevolume / 1000);  
%each step is a constant fraction of the overall MO volume
mass_liquid     = Mantlemass-mass_solidified;
mass_solidified = mass_solidified + mass_this_shell;
new_mass_liquid = Mantlemass-mass_solidified;


if P(j) < Layer6P;
    new_liq = liquid(j-1,:);
else        % Note that solid composition includes interstitial liquids
    %    new_liq = liquid(j-1,:)*(liqfrac(j)/liqfrac(j-1)) - solid(j,:)*(1-liqfrac(j)/liqfrac(j-1))*(1-intliqx) - liquid(j-1,:)*(1-liqfrac(j)/liqfrac(j-1))*(intliqx);
    %   new_liq = liquid(j-1,:) - solid(j,:)*(1-liqfrac(j)/liqfrac(j-1)); %*(1-intliqx); % - liquid(j-1,:)*(1-liqfrac(j)/liqfrac(j-1))*(intliqx);
    %new_liq = liquid(j-1,:) - solid(j,:)*((Mantlemass*0.001)/(liqfrac(j)*Mantlemass)); % - liquid(j-1,:)*(intliqx)*((Mantlemass*0.001)/(liqfrac(j)*Mantlemass));
    new_liq = (liquid(j-1,:) * mass_liquid - mass_this_shell * solid(j,:))/new_mass_liquid;
    if min(new_liq) < 0
        disp('new_liq < 0')
    end
    liq_comp = 100*new_liq./(sum(new_liq));  % renormalized liquid composition
end



liquid(j,:) = liq_comp;
MG(j,1) = (liq_comp(4)/40.311)./(liq_comp(4)/40.311 + liq_comp(3)/71.846);
MG(j,2) = (solid(j,4)/40.311)./(solid(j,4)/40.311 + solid(j,3)/71.846);
D(j) = density;
Do(j) = densityzero;
Dsol(j) = densitysol;

if liquid(j,1) <= 3;
    liquid(j,1) = 3;
    disp(['No Si at step: ', num2str(j)]);
    liquid(j,:) = 100*liquid(j,:)/sum(liquid(j,:));
end
if liquid(j,2) <= 2;
    liquid(j,1) = 2;
    disp(['No Al at step: ', num2str(j)]);
    liquid(j,:) = 100*liquid(j,:)/sum(liquid(j,:));
end
if liquid(j,3) <= 1; liquid(j,3) = 1;
    disp(['No Fe at step: ', num2str(j)]);
    liquid(j,:) = 100*liquid(j,:)/sum(liquid(j,:));
end
if liquid(j,4) <= 1;
    liquid(j,4) = 1;
    disp(['No Mg at step: ', num2str(j)]);
    liquid(j,:) = 100*liquid(j,:)/sum(liquid(j,:)); % normalization
end
if liquid(j,5) <= 1;
    liquid(j,5) = 1;
    disp(['No Ca at step: ', num2str(j)]);
    liquid(j,:) = 100*liquid(j,:)/sum(liquid(j,:));
end

%        if liquid(j-1,6) <= 1e-7; liquid(j-1,5) = 1e-6; disp(['No Sm at step: ', num2str(j)]); liquid(j,:) = 100*liquid(j,:)/sum(liquid(j,:));end
%        if liquid(j-1,7) <= 1e-7; liquid(j-1,5) = 1e-6; disp(['No Nd at step: ', num2str(j)]); liquid(j,:) = 100*liquid(j,:)/sum(liquid(j,:));end
%        if liquid(j-1,8) <= 1e-7; liquid(j-1,5) = 1e-6; disp(['No Lu at step: ', num2str(j)]); liquid(j,:) = 100*liquid(j,:)/sum(liquid(j,:));end
%        if liquid(j-1,9) <= 1e-7; liquid(j-1,5) = 1e-6; disp(['No Hf at step: ', num2str(j)]); liquid(j,:) = 100*liquid(j,:)/sum(liquid(j,:));end
%        if P(j)<= 1       % do not evolve liquid composition any more - can't handle it
%             liquid(j,1:9) = liquid(j-1,1:9); liquid(j,:) = 100*liquid(j,:)/sum(liquid(j,:));
%        end

