% fractionatedeepEarth.m

% in km and GPa
% in solid and liquid matrices:
%      1 = SiO2, 2 = Al2O3, 3 = FeO, 4 = MgO, 5 = CaO, 6 = Sm, 7 = Nd, 8 = Lu, 9 = Hf, 10 = OH, 11 = C
% VECTORS: D = density at r; Do = density at 1 atm and 298 K; M = mass interior to r, P = pressure at r; MG = Mg# liquid, solid

%change radius into pressure GPa

P(j) = -0.0374*(r(j)/1000) + 238.5372;   %***[GPa] from radius in km, from 130GPa at 2900 depth (CMB) to 0 at 6378 (surface)

%%
% %%%%%% LAYER 0 %%%%%%%%%%%%%%%%%%%% post-perovskite
%display('Layer 0: Post-perovskite')
if (P(j) > Layer0P);
   intliqx = intliq0;
   [eqminppv, Mgnumppv] = postperovskite(liq_comp);
   
   ppvdensity = postperovskitedensity(Mgnumppv, 0, P(j), Tsolid(j));
   ppvdensityzero = postperovskitedensity(Mgnumppv, 0, 1e-4, 1);
   ppvdensitysol = postperovskitedensity(Mgnumppv, 0, 1e-4, Tsolid(j));
   
   
   intliq(j,:) = liquid(j-1,:);
   [intliqdensity, intliqudensityzero, intliqdensitysol] = NewBirchMurnliquid(intliq(j,:), Tsolid(j), P(j));
   
   solid(j,:) = ppv0*eqminppv + intliq0*intliq(j,:);
   solidOH(j,:) = ppv0*eqminppv(10);
   intliqOH(j,:) = intliq0*intliq(j,10);
   
   Mgratio = ((liq_comp(4)/40.311)./(liq_comp(4)/40.311 + liq_comp(3)/71.846))/((solid(j,4)/40.311)./(solid(j,4)/40.311 + solid(j,3)/71.846));
   density = ppv0*ppvdensity + intliq0*Mgratio*ppvdensity;
   densityzero = ppv0*ppvdensityzero + intliq0*Mgratio*ppvdensityzero;
   densitysol = ppv0*ppvdensitysol + intliq0*Mgratio*ppvdensitysol;
   
   %%%%% ????
   pperov(j) = ppvdensity; magnesio(j) = 0; gam(j) = 0; bet(j) = 0;
   maj(j) = 0; stish(j) = 0; cpyrox(j) = 0; opyrox(j) = 0;
   gar(j) = 0; oliv(j) = 0; iliq(j) = intliqdensity;
   
   finallayercalculations
   marker0 = j-1;
   
   return
   
end

%%
% %%%%%% LAYER 1 %%%%%%%%%%%%%%%%%%%% mw + perov
if P(j) > Layer1P;
%display('Layer 1: Perovskite, Magnesiowustite')
    intliqx = intliq1;
    %    disp(['Layer 2, index = ', num2str(j)])
    [eqminper, Mgnumper, Perc_Al, Perc_Ca, Perc_MgFe] = perovskite(liq_comp);
    pdensity = perovskitedensity(Mgnumper, Perc_Al, Perc_Ca, Perc_MgFe, P(j), Tsolid(j));
    pdensityzero = perovskitedensity(Mgnumper, Perc_Al, Perc_Ca, Perc_MgFe, 1e-4, 1);
    pdensitysol = perovskitedensity(Mgnumper, Perc_Al, Perc_Ca, Perc_MgFe, 1e-4, Tsolid(j));
    
    [Mgnummw, eqminmw] = magnesiowustite(liq_comp);
    mwdensity = magnesiowustitedensity(Mgnummw, P(j), Tsolid(j));
    mwdensityzero = magnesiowustitedensity(Mgnummw, 1e-4, 1);
    mwdensitysol = magnesiowustitedensity(Mgnummw, 1e-4, Tsolid(j));

    intliq(j,:) = liquid(j-1,:);
    [intliqdensity, intliqdensityzero, intliqdensitysol]=NewBirchMurnliquid(intliq(j,:), Tsolid(j), P(j));
    %density of solidified interstitial liquid is assumed to be related to the rest of the solids by the ratio of their MG#s,
    %and ignores the fact it would crystallize other phases
    solid(j,:) = per1*eqminper + mw1*eqminmw + intliq1*intliq(j,:);
    solidOH(j,:) = per1*eqminper(10) + mw1*eqminmw(10);
    intliqOH(j,:) = intliq1*intliq(j,10);
    Mgratio = [(liq_comp(4)/40.311)./(liq_comp(4)/40.311 + liq_comp(3)/71.846)]/[(solid(j,4)/40.311)./(solid(j,4)/40.311 + solid(j,3)/71.846)];
    density = per1*pdensity + mw1*mwdensity + intliq1*Mgratio*[(per1/(per1 + mw1))*pdensity + (mw1/(per1 + mw1))*mwdensity];
    densityzero = per1*pdensityzero + mw1*mwdensityzero + intliq1*Mgratio*[(per1/(per1 + mw1))*pdensityzero + (mw1/(per1 + mw1))*mwdensityzero];
    densitysol = per1*pdensitysol + mw1*mwdensitysol + intliq1*Mgratio*[(per1/(per1 + mw1))*pdensitysol + (mw1/(per1 + mw1))*mwdensitysol];

    perov(j) = pdensity; magnesio(j) = mwdensity; gam(j) = 0; bet(j) = 0;
    maj(j) = 0; stish(j) = 0; cpyrox(j) = 0; opyrox(j) = 0;
    gar(j) = 0; oliv(j) = 0; iliq(j) = intliqdensity;

    finallayercalculations
    marker1 = j-1;
    return
end

% %%%%%% LAYER 2 %%%%%%%%%%%%%%%%%%%% gammaspinel + majorite
if P(j) > Layer2P;
%display('Layer 2: Gamma spinel, Majorite')
    intliqx = intliq2;
    %    disp(['Layer 2, index = ', num2str(j)])
    [eqmingamma, Mgnumg] = gamma(liq_comp);
    gdensity = gammadensity(Mgnumg, P(j), Tsolid(j));
    gdensityzero = gammadensity(Mgnumg, 1e-4, 1);
    gdensitysol = gammadensity(Mgnumg, 1e-4, Tsolid(j));

    [eqminmaj, Mgnummaj, Perc_Ca, Perc_AlMgFe, Perc_MgFe] = majorite(liq_comp);
    majdensity = majoritedensity(Mgnummaj, Perc_Ca, Perc_AlMgFe, Perc_MgFe, P(j), Tsolid(j));
    majdensityzero = majoritedensity(Mgnummaj, Perc_Ca, Perc_AlMgFe, Perc_MgFe, 1e-4, 1);
    majdensitysol = majoritedensity(Mgnummaj, Perc_Ca, Perc_AlMgFe, Perc_MgFe, 1e-4, Tsolid(j));

    intliq(j,:) = liquid(j-1,:);
    [intliqdensity, intliqdensityzero, intliqdensitysol]=NewBirchMurnliquid(intliq(j,:), Tsolid(j), P(j));
    %density of solidified interstitial liquid is assumed to be related to the rest of the solids by the ratio of their MG#s,
    %and ignores the fact it would crystallize other phases
    solid(j,:) = gamma2*eqmingamma + maj2*eqminmaj + intliq2*intliq(j,:);
    solidOH(j,:) = gamma2*eqmingamma(10) + maj2*eqminmaj(10);
    intliqOH(j,:) = intliq2*intliq(j,10);
    Mgratio = [(liq_comp(4)/40.311)./(liq_comp(4)/40.311 + liq_comp(3)/71.846)]/[(solid(j,4)/40.311)./(solid(j,4)/40.311 + solid(j,3)/71.846)];
    density = gamma2*gdensity + maj2*majdensity + intliq2*Mgratio*[(gamma2/(gamma2 + maj2))*gdensity + (maj2/(gamma2 + maj2))*majdensity];
    densityzero = gamma2*gdensityzero + maj2*majdensityzero + intliq2*Mgratio*[(gamma2/(gamma2 + maj2))*gdensityzero + (maj2/(gamma2 + maj2))*majdensityzero];
    densitysol = gamma2*gdensitysol + maj2*majdensitysol + intliq2*Mgratio*[(gamma2/(gamma2 + maj2))*gdensitysol + (maj2/(gamma2 + maj2))*majdensitysol];

    perov(j) = 0; magnesio(j) = 0; gam(j) = gdensity;  bet(j) = 0;
    maj(j) = majdensity; stish(j) = 0; cpyrox(j) = 0; opyrox(j) = 0;
    gar(j) = 0; oliv(j) = 0; iliq(j) = intliqdensity;

    finallayercalculations
    marker2 = j-1;
    return
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LAYER 3 %%%%%%%%%%%%%%%%%%%% gammaspinel + majorite
if P(j) > Layer3P;
%display('Layer 3: Gamma spinel, Majorite')
    intliqx = intliq3;
    %    disp(['Layer 3, index = ', num2str(j)])
    [eqminbeta, Mgnumb] = beta(liq_comp);
    bdensity = betadensity(Mgnumb, P(j), Tsolid(j));
    bdensityzero = betadensity(Mgnumb, 1e-4, 1);
    bdensitysol = betadensity(Mgnumb, 1e-4, Tsolid(j));

    [eqminmaj, Mgnummaj, Perc_Ca, Perc_AlMgFe, Perc_MgFe] = majorite(liq_comp);
    majdensity = majoritedensity(Mgnummaj, Perc_Ca, Perc_AlMgFe, Perc_MgFe, P(j), Tsolid(j));
    majdensityzero = majoritedensity(Mgnummaj, Perc_Ca, Perc_AlMgFe, Perc_MgFe, 1e-4, 1);
    majdensitysol = majoritedensity(Mgnummaj, Perc_Ca, Perc_AlMgFe, Perc_MgFe, 1e-4, Tsolid(j));

    [eqmincpx, Mgnumcpx, CaMg] = clinopyroxene(liq_comp);
    cpxdensity = clinopyroxenedensity(Mgnumcpx, CaMg, P(j), Tsolid(j));
    cpxdensityzero = clinopyroxenedensity(Mgnumcpx, CaMg, 1e-4, 1);
    cpxdensitysol = clinopyroxenedensity(Mgnumcpx, CaMg, 1e-4, Tsolid(j));

    intliq(j,:) = liquid(j-1,:);
    [intliqdensity, intliqdensityzero, intliqdensitysol]=NewBirchMurnliquid(intliq(j,:), Tsolid(j), P(j));
    %density of solidified interstitial liquid is assumed to be related to the rest of the solids by the ratio of their MG#s,
    %and ignores the fact it would crystallize other phases
    solid(j,:) = beta3*eqminbeta + maj3*eqminmaj + cpx3*eqmincpx + intliq3*intliq(j,:);
    solidOH(j,:) = beta3*eqminbeta(10) + maj3*eqminmaj(10) + cpx3*eqmincpx(10);
    intliqOH(j,:) = intliq3*intliq(j,10);
    Mgratio = [(liq_comp(4)/40.311)./(liq_comp(4)/40.311 + liq_comp(3)/71.846)]/[(solid(j,4)/40.311)./(solid(j,4)/40.311 + solid(j,3)/71.846)];
    density = beta3*bdensity + maj3*majdensity + cpx3*cpxdensity + intliq3*Mgratio*[(beta3/(beta3 + maj3 + cpx3))*bdensity + (maj3/(beta3 + maj3 + cpx3))*majdensity + (cpx3/(beta3 + maj3 + cpx3))*cpxdensity];
    densityzero = beta3*bdensityzero + maj3*majdensityzero + cpx3*cpxdensityzero + intliq3*Mgratio*[(beta3/(beta3 + maj3 + cpx3))*bdensityzero + (maj3/(beta3 + maj3 + cpx3))*majdensityzero + (cpx3/(beta3 + maj3 + cpx3))*cpxdensityzero];
    densitysol = beta3*bdensitysol + maj3*majdensitysol + cpx3*cpxdensitysol + intliq3*Mgratio*[(beta3/(beta3 + maj3 + cpx3))*bdensitysol + (maj3/(beta3 + maj3 + cpx3))*majdensitysol + (cpx3/(beta3 + maj3 + cpx3))*cpxdensitysol];

    perov(j) = 0; magnesio(j) = 0; bet(j) = bdensity; gam(j) = 0;
    maj(j) = majdensity; stish(j) = 0; cpyrox(j) = cpxdensity; opyrox(j) = 0;
    gar(j) = 0; oliv(j) = 0; iliq(j) = intliqdensity;

    finallayercalculations

    return
end

%
% % % %%%%%% LAYER 4 BELOW SEPTUM %%%%%%%%%%%%%%%%%%%% garnet only
% if P(j) > Layer4PA;
%   intliqx = intliq4;
% if liquid(j-1,2) > 1.5;	% only garnet falls out and is sequestered; ol + pyx rehomogenize and are therefore not calculated here
% % disp(['Layer 4, garnet fractionation, index = ', num2str(j)])
%        [eqmingar, Mgnumgar, Perc_AlCa, Perc_AlMgFe] = garnet(liq_comp);
%        gardensity = garnetdensity(Mgnumgar, Perc_AlCa, Perc_AlMgFe, P(j), Tsolid(j));
%        gardensityzero = garnetdensity(Mgnumgar, Perc_AlCa, Perc_AlMgFe, 1e-4, 1);
%        gardensitysol = garnetdensity(Mgnumgar, Perc_AlCa, Perc_AlMgFe, 1e-4, Tsolid(j));
%
%        intliq(j,:) = liquid(j-1,:);
%        [intliqdensity, intliqdensityzero, intliqdensitysol]=NewBirchMurnliquid(intliq(j,:), Tsolid(j), P(j));
%        %density of solidified interstitial liquid is assumed to be related to the rest of the solids by the ratio of their MG#s,
%             %and ignores the fact it would crystallize other phases
%        solid(j,:) = (1-intliq4)*eqmingar + intliq4*intliq(j,:);
%        Mgratio = [(liq_comp(4)/40.311)./(liq_comp(4)/40.311 + liq_comp(3)/71.846)]/[(solid(j,4)/40.311)./(solid(j,4)/40.311 + solid(j,3)/71.846)];
%        density = (1-intliq4)*gardensity + intliq4*Mgratio*gardensity;
%        densityzero = (1-intliq4)*gardensityzero + intliq4*Mgratio*gardensityzero;
%        densitysol = (1-intliq4)*gardensitysol + intliq4*Mgratio*gardensitysol;
%
%        perov(j) = 0; magnesio(j) = 0; gam(j) = 0;  bet(j) = 0;
%        maj(j) = 0; stish(j) = 0; pyrox(j) = 0;
%        gar(j) = gardensity; oliv(j) = 0; iliq(j) = intliqdensity;
%
% finallayercalculations
%
% if liquid(j-1,2) < 2; Layer4PA = P(j); end  % will not go back to garnet once Al exhausted
%
% return
% end
% end
%
% %disp('End of garnet fractionation - no more Al'); j

%% %%%%%% LAYER 4 %%%%%%%%%%%%%%%%%%%% alpha + opx + cpx + gt
if P(j) > Layer4P;
%display('Layer 4: Garnet, Alpha olivine, Clinopyroxene, Orthopyroxene')
    intliqx = intliq4;
    %    disp(['Layer 4 after garnet, index = ', num2str(j)])
    [eqmingar, Mgnumgar, Perc_AlCa, Perc_AlMgFe] = garnet(liq_comp);
    gardensity = garnetdensity(Mgnumgar, Perc_AlCa, Perc_AlMgFe, P(j), Tsolid(j));
    gardensityzero = garnetdensity(Mgnumgar, Perc_AlCa, Perc_AlMgFe, 1e-4, 1);
    gardensitysol = garnetdensity(Mgnumgar, Perc_AlCa, Perc_AlMgFe, 1e-4, Tsolid(j));

    [eqminalpha, Mgnumalpha] = olivine(liq_comp);
    alphadensity = olivinedensity(Mgnumalpha, P(j), Tsolid(j));
    alphadensityzero = olivinedensity(Mgnumalpha, 1e-4, 1);
    alphadensitysol = olivinedensity(Mgnumalpha, 1e-4, Tsolid(j));

    [eqmincpx, Mgnumcpx, CaMg] = clinopyroxene(liq_comp);
    cpxdensity = clinopyroxenedensity(Mgnumcpx, CaMg, P(j), Tsolid(j));
    cpxdensityzero = clinopyroxenedensity(Mgnumcpx, CaMg, 1e-4, 1);
    cpxdensitysol = clinopyroxenedensity(Mgnumcpx, CaMg, 1e-4, Tsolid(j));

    [eqminopx, Mgnumopx, CaMg] = orthopyroxene(liq_comp);
    opxdensity = orthopyroxenedensity(Mgnumopx, CaMg, P(j), Tsolid(j));
    opxdensityzero = orthopyroxenedensity(Mgnumopx, CaMg, 1e-4, 1);
    opxdensitysol = orthopyroxenedensity(Mgnumopx, CaMg, 1e-4, Tsolid(j));

    intliq(j,:) = liquid(j-1,:);
    [intliqdensity, intliqdensityzero, intliqdensitysol]=NewBirchMurnliquid(intliq(j,:), Tsolid(j), P(j));
    %density of solidified interstitial liquid is assumed to be related to the rest of the solids by the ratio of their MG#s,
    %and ignores the fact it would crystallize other phases
    solid(j,:) = gar4*eqmingar + alpha4*eqminalpha + cpx4*eqmincpx + opx4*eqminopx + intliq4*intliq(j,:);
    solidOH(j,:) = gar4*eqmingar(10) + alpha4*eqminalpha(10) + cpx4*eqmincpx(10) + opx4*eqminopx(10);
    intliqOH(j,:) = intliq4*intliq(j,10);
    Mgratio = [(liq_comp(4)/40.311)./(liq_comp(4)/40.311 + liq_comp(3)/71.846)]/[(solid(j,4)/40.311)./(solid(j,4)/40.311 + solid(j,3)/71.846)];
    density = gar4*gardensity + alpha4*alphadensity + cpx4*cpxdensity + opx4*opxdensity + intliq4*Mgratio*[(alpha4/(alpha4 + cpx4))*alphadensity + (cpx4/(alpha4 + cpx4))*cpxdensity];
    densityzero = gar4*gardensityzero + alpha4*alphadensityzero + cpx4*cpxdensityzero + opx4*opxdensityzero + intliq4*Mgratio*[(alpha4/(alpha4 + cpx4))*alphadensityzero + (cpx4/(alpha4 + cpx4))*cpxdensityzero];
    densitysol = gar4*gardensitysol + alpha4*alphadensitysol + cpx4*cpxdensitysol + opx4*opxdensitysol + intliq4*Mgratio*[(alpha4/(alpha4 + cpx4))*alphadensitysol + (cpx4/(alpha4 + cpx4))*cpxdensitysol];

    perov(j) = 0; magnesio(j) = 0; bet(j) = 0; gam(j) = 0;  spi(j) = 0;  plag(j) = 0;
    maj(j) = 0; stish(j) = 0; cpyrox(j) = cpxdensity; opyrox(j) = opxdensity;
    gar(j) = gardensity; oliv(j) = alphadensity; iliq(j) = intliqdensity;

    finallayercalculations

    return
end

%%
% % %%%%%% LAYER 5 %%%%%%%%%%%%%%%%%%%% spin + opx + cpx + alpha
if P(j) > Layer5P;
%display('Layer 5: Spinel, Alpha olivine, Clinopyroxene, Orthopyroxene')
    intliqx = intliq5;
    %    disp(['Layer 5, index = ', num2str(j)])
    [eqminspin, Mgnumspin] = spinel(liq_comp);
    spindensity = spineldensity(Mgnumspin, P(j), Tsolid(j));
    spindensityzero = spineldensity(Mgnumspin, 1e-4, 1);
    spindensitysol = spineldensity(Mgnumspin, 1e-4, Tsolid(j));

    [eqminalpha, Mgnumalpha] = olivine(liq_comp);
    alphadensity = olivinedensity(Mgnumalpha, P(j), Tsolid(j));
    alphadensityzero = olivinedensity(Mgnumalpha, 1e-4, 1);
    alphadensitysol = olivinedensity(Mgnumalpha, 1e-4, Tsolid(j));

    [eqmincpx, Mgnumcpx, CaMg] = clinopyroxene(liq_comp);
    cpxdensity = clinopyroxenedensity(Mgnumcpx, CaMg, P(j), Tsolid(j));
    cpxdensityzero = clinopyroxenedensity(Mgnumcpx, CaMg, 1e-4, 1);
    cpxdensitysol = clinopyroxenedensity(Mgnumcpx, CaMg, 1e-4, Tsolid(j));

    [eqminopx, Mgnumopx, CaMg] = orthopyroxene(liq_comp);
    opxdensity = orthopyroxenedensity(Mgnumopx, CaMg, P(j), Tsolid(j));
    opxdensityzero = orthopyroxenedensity(Mgnumopx, CaMg, 1e-4, 1);
    opxdensitysol = orthopyroxenedensity(Mgnumopx, CaMg, 1e-4, Tsolid(j));

    intliq(j,:) = liquid(j-1,:);
    [intliqdensity, intliqdensityzero, intliqdensitysol]=NewBirchMurnliquid(intliq(j,:), Tsolid(j), P(j));
    %density of solidified interstitial liquid is assumed to be related to the rest of the solids by the ratio of their MG#s,
    %and ignores the fact it would crystallize other phases
    solid(j,:) = spin5*eqminspin + alpha5*eqminalpha + cpx5*eqmincpx + opx5*eqminopx + intliq5*intliq(j,:);
    solidOH(j,:) = spin5*eqminspin(10) + alpha5*eqminalpha(10) + cpx5*eqmincpx(10) + opx5*eqminopx(10);
    intliqOH(j,:) = intliq5*intliq(j,10);
    Mgratio = [(liq_comp(4)/40.311)./(liq_comp(4)/40.311 + liq_comp(3)/71.846)]/[(solid(j,4)/40.311)./(solid(j,4)/40.311 + solid(j,3)/71.846)];
    density = spin5*spindensity + alpha5*alphadensity + cpx5*cpxdensity + opx5*opxdensity + intliq5*Mgratio*[(alpha5/(alpha5 + cpx5))*alphadensity + (cpx5/(alpha5 + cpx5))*cpxdensity];
    densityzero = spin5*spindensityzero + alpha5*alphadensityzero + cpx5*cpxdensityzero + opx5*opxdensityzero + intliq5*Mgratio*[(alpha5/(alpha5 + cpx5))*alphadensityzero + (cpx5/(alpha5 + cpx5))*cpxdensityzero];
    densitysol = spin5*spindensitysol + alpha5*alphadensitysol + cpx5*cpxdensitysol + opx5*opxdensitysol + intliq5*Mgratio*[(alpha5/(alpha5 + cpx5))*alphadensitysol + (cpx5/(alpha5 + cpx5))*cpxdensitysol];

    perov(j) = 0; magnesio(j) = 0; bet(j) = 0; gam(j) = 0; spi(j) = spindensity;  plag(j) = 0;
    maj(j) = 0; stish(j) = 0; cpyrox(j) = cpxdensity; opyrox(j) = opxdensity;
    gar(j) = 0; oliv(j) = alphadensity; iliq(j) = intliqdensity;

    finallayercalculations

    return
end

%%
%%%%%%%%%%%%%%%%%% LAYER 6 %%%%%%%%%%%%%%%%%%%% plag + opx + cpx + alpha
if P(j) > Layer6P;
%display('Layer 6: Plagioclase, Alpha olivine, Clinopyroxene, Orthopyroxene')
    intliqx = intliq6;
    %    disp(['Layer 6, index = ', num2str(j)])
    [eqminplag, Anplag] = plagioclase(liq_comp);
    plagdensity = plagioclasedensity(Anplag, P(j), Tsolid(j));
    plagdensityzero = plagioclasedensity(Anplag, 1e-4, 1);
    plagdensitysol = plagioclasedensity(Anplag, 1e-4, Tsolid(j));

    [eqminalpha, Mgnumalpha] = olivine(liq_comp);
    alphadensity = olivinedensity(Mgnumalpha, P(j), Tsolid(j));
    alphadensityzero = olivinedensity(Mgnumalpha, 1e-4, 1);
    alphadensitysol = olivinedensity(Mgnumalpha, 1e-4, Tsolid(j));

    [eqmincpx, Mgnumpyx, CaMg] = clinopyroxene(liq_comp);
    cpxdensity = clinopyroxenedensity(Mgnumpyx, CaMg, P(j), Tsolid(j));
    cpxdensityzero = clinopyroxenedensity(Mgnumpyx, CaMg, 1e-4, 1);
    cpxdensitysol = clinopyroxenedensity(Mgnumpyx, CaMg, 1e-4, Tsolid(j));

    [eqminopx, Mgnumpyx, CaMg] = orthopyroxene(liq_comp);
    opxdensity = orthopyroxenedensity(Mgnumpyx, CaMg, P(j), Tsolid(j));
    opxdensityzero = orthopyroxenedensity(Mgnumpyx, CaMg, 1e-4, 1);
    opxdensitysol = orthopyroxenedensity(Mgnumpyx, CaMg, 1e-4, Tsolid(j));

    intliq(j,:) = liquid(j-1,:);
    %        [intliqdensity, intliqdensityzero, intliqdensitysol]=NewBirchMurnliquid(intliq(j,:), Tsolid(j), P(j));
    intliqdensitysol = 1000*surfacemeltdensity(intliq(j,:), Tsolid(j)); % Birch-Murnaghan solution breaks down at low P
    intliqdensity = intliqdensitysol+100*P(j);
    intliqdensityzero = 1000*surfacemeltdensity(intliq(j,:), 1);
    %density of solidified interstitial liquid is assumed to be related to the rest of the solids by the ratio of their MG#s,
    %and ignores the fact it would crystallize other phases
    solid(j,:) = plag6*eqminplag + alpha6*eqminalpha + cpx6*eqmincpx + opx6*eqminopx + intliq6*intliq(j,:);
    solidOH(j,:) = plag6*eqminplag(10) + alpha6*eqminalpha(10) + cpx6*eqmincpx(10) + opx6*eqminopx(10);
    intliqOH(j,:) = intliq6*intliq(j,10);
    Mgratio = [(liq_comp(4)/40.311)./(liq_comp(4)/40.311 + liq_comp(3)/71.846)]/[(solid(j,4)/40.311)./(solid(j,4)/40.311 + solid(j,3)/71.846)];
    density = plag6*plagdensity + alpha6*alphadensity + cpx6*cpxdensity + opx6*opxdensity + intliq6*Mgratio*[(alpha6/(alpha6 + cpx6))*alphadensity + (cpx6/(alpha6 + cpx6))*cpxdensity];
    densityzero = plag6*plagdensityzero + alpha6*alphadensityzero + cpx6*cpxdensityzero + opx6*opxdensityzero + intliq6*Mgratio*[(alpha6/(alpha6 + cpx6))*alphadensityzero + (cpx6/(alpha6 + cpx6))*cpxdensityzero];
    densitysol = plag6*plagdensitysol + alpha6*alphadensitysol + cpx6*cpxdensitysol + opx6*opxdensitysol + intliq6*Mgratio*[(alpha6/(alpha6 + cpx6))*alphadensitysol + (cpx6/(alpha6 + cpx6))*cpxdensitysol];

    perov(j) = 0; magnesio(j) = 0; bet(j) = 0; gam(j) = 0; spi(j) = 0; plag(j) = plagdensity;
    maj(j) = 0; stish(j) = 0; cpyrox(j) = cpxdensity; opyrox(j) = opxdensity;
    gar(j) = 0; oliv(j) = alphadensity; iliq(j) = intliqdensity;

    finallayercalculations

    return
end

%%
%%%%%%%%%%%%%%%%%% LAYER 7 %%%%%%%%%%%% same as layer 6; fixed liquid
if P(j) > Layer7P;
%display('Layer 7: Same as layer 6 (fixed liquid)')
    intliqx = intliq7;
    %    disp(['Layer 7, index = ', num2str(j)])
    [eqminplag, Anplag] = plagioclase(liq_comp);
    plagdensity = plagioclasedensity(Anplag, P(j), Tsolid(j));
    plagdensityzero = plagioclasedensity(Anplag, 1e-4, 1);
    plagdensitysol = plagioclasedensity(Anplag, 1e-4, Tsolid(j));

    [eqminalpha, Mgnumalpha] = olivine(liq_comp);
    alphadensity = olivinedensity(Mgnumalpha, P(j), Tsolid(j));
    alphadensityzero = olivinedensity(Mgnumalpha, 1e-4, 1);
    alphadensitysol = olivinedensity(Mgnumalpha, 1e-4, Tsolid(j));

    [eqmincpx, Mgnumpyx, CaMg] = clinopyroxene(liq_comp);
    cpxdensity = clinopyroxenedensity(Mgnumpyx, CaMg, P(j), Tsolid(j));
    cpxdensityzero = clinopyroxenedensity(Mgnumpyx, CaMg, 1e-4, 1);
    cpxdensitysol = clinopyroxenedensity(Mgnumpyx, CaMg, 1e-4, Tsolid(j));

    [eqminopx, Mgnumpyx, CaMg] = orthopyroxene(liq_comp);
    opxdensity = orthopyroxenedensity(Mgnumpyx, CaMg, P(j), Tsolid(j));
    opxdensityzero = orthopyroxenedensity(Mgnumpyx, CaMg, 1e-4, 1);
    opxdensitysol = orthopyroxenedensity(Mgnumpyx, CaMg, 1e-4, Tsolid(j));

    intliq(j,:) = liquid(j-1,:);
    %        [intliqdensity, intliqdensityzero, intliqdensitysol]=NewBirchMurnliquid(intliq(j,:), Tsolid(j), P(j));
    intliqdensitysol = 1000*surfacemeltdensity(intliq(j,:), Tsolid(j)); % Birch-Murnaghan solution breaks down at low P
    intliqdensity = intliqdensitysol+100*P(j);
    intliqdensityzero = 1000*surfacemeltdensity(intliq(j,:), 1);
    %density of solidified interstitial liquid is assumed to be related to the rest of the solids by the ratio of their MG#s,
    %and ignores the fact it would crystallize other phases
    solid(j,:) = plag7*eqminplag + alpha7*eqminalpha + cpx7*eqmincpx + opx7*eqminopx + intliq7*intliq(j,:);
    solidOH(j,:) = plag7*eqminplag(10) + alpha7*eqminalpha(10) + cpx7*eqmincpx(10) + opx7*eqminopx(10);
    intliqOH(j,:) = intliq7*intliq(j,10);
    Mgratio = [(liq_comp(4)/40.311)./(liq_comp(4)/40.311 + liq_comp(3)/71.846)]/[(solid(j,4)/40.311)./(solid(j,4)/40.311 + solid(j,3)/71.846)];
    density = plag7*plagdensity + alpha7*alphadensity + cpx7*cpxdensity + opx7*opxdensity + intliq7*Mgratio*[(alpha7/(alpha7 + cpx7))*alphadensity + (cpx7/(alpha7 + cpx7))*cpxdensity];
    densityzero = plag7*plagdensityzero + alpha7*alphadensityzero + cpx7*cpxdensityzero + opx7*opxdensityzero + intliq7*Mgratio*[(alpha7/(alpha7 + cpx7))*alphadensityzero + (cpx7/(alpha7 + cpx7))*cpxdensityzero];
    densitysol = plag7*plagdensitysol + alpha7*alphadensitysol + cpx7*cpxdensitysol + opx7*opxdensitysol + intliq7*Mgratio*[(alpha7/(alpha7 + cpx7))*alphadensitysol + (cpx7/(alpha7 + cpx7))*cpxdensitysol];

    perov(j) = 0; magnesio(j) = 0; bet(j) = 0; gam(j) = 0; spi(j) = 0; plag(j) = plagdensity;
    maj(j) = 0; stish(j) = 0; cpyrox(j) = cpxdensity; opyrox(j) = opxdensity;
    gar(j) = 0; oliv(j) = alphadensity; iliq(j) = intliqdensity;

    finallayercalculations

    return
end