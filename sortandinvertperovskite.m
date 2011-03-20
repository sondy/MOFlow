% sortandinvertperovskite.m

disp('In sortandinvertperovskite')

% includes separate recalculation of density for perovskite + mw layer 5/9/2007
% altered for integrated code on 09/23/2006
% takes results of any profile and inverts it by sorting by density
% first calculates density of first layer as mineralogy of second layer
% plots new mineralogy densities in figures 1 and 2
% for new planet, change two lines of "solidus" below

endi = j;

% for i = 1:endi      % make vectors of mass of mantle only, and of new radii
%     MantleM(i) = M(i) - M(1);
% end
%
% Mean = mean(MantleM);
% MiddleMass = find(MantleM >= Mean);   % makes marker at index of half of mantle mass
%
%

%%%%%this loop recalculates perovskite + magnesiowustite as olivine as cpx + olivine + garnet
if exist('marker1','var')  %if the ocean is shallow, marker1  and marker2 may not get defined.

    for i = 1:marker1		% don't actually calculate the new compositions, but
        % estimate densities based on usual parameters
        Mgnumpyx = MG(i,2) + 0.01;
        CaMg = 0.2;
        cpxdensity = clinopyroxenedensity(Mgnumpyx, CaMg, P(i), Tsolid(i));
        cpxdensityzero = clinopyroxenedensity(Mgnumpyx, CaMg, 1e-4, 1);
        cpxdensitysol = clinopyroxenedensity(Mgnumpyx, CaMg, 1e-4, Tsolid(i));

        Mgnumgar = MG(i,2) - 0.02;
        Perc_AlCa = 0.05;
        Perc_AlMgFe = 0.95;
        gardensity = garnetdensity(Mgnumgar, Perc_AlCa, Perc_AlMgFe, P(i), Tsolid(i));
        gardensityzero = garnetdensity(Mgnumgar, Perc_AlCa, Perc_AlMgFe, 1e-4, 1);
        gardensitysol = garnetdensity(Mgnumgar, Perc_AlCa, Perc_AlMgFe, 1e-4, Tsolid(i));

        Mgnumalpha = MG(i,2);
        alphadensity = olivinedensity(Mgnumalpha, P(i), Tsolid(i));
        aldensityzero = olivinedensity(Mgnumalpha, 1e-4, 1);
        aldensitysol = olivinedensity(Mgnumalpha, 1e-4, Tsolid(i));

        density = 0.04*cpxdensity + 0.0225*gardensity + 0.935*alphadensity;
        densityzero = 0.04*cpxdensityzero + 0.0225*gardensityzero + 0.935*aldensityzero;
        densitysol = 0.04*cpxdensitysol + 0.0225*gardensitysol + 0.935*aldensitysol;

        D(i) = density;
        Do(i) = densityzero;
        Dsol(i) = densitysol;
    end
else
    marker1 = 1;
end

if exist('marker2','var')  %if the ocean is shallow, marker1  and marker2 may not get defined.
    %%%%%this loop recalculates majorite + gamma olivine as cpx + olivine + garnet
    for i = marker1:marker2		% don't actually calculate the new compositions, but
        % estimate densities based on usual parameters
        Mgnumpyx = MG(i,2) + 0.01;
        CaMg = 0.2;
        cpxdensity = clinopyroxenedensity(Mgnumpyx, CaMg, P(i), Tsolid(i));
        cpxdensityzero = clinopyroxenedensity(Mgnumpyx, CaMg, 1e-4, 1);
        cpxdensitysol = clinopyroxenedensity(Mgnumpyx, CaMg, 1e-4, Tsolid(i));

        Mgnumgar = MG(i,2) - 0.02;
        Perc_AlCa = 0.05;
        Perc_AlMgFe = 0.95;
        gardensity = garnetdensity(Mgnumgar, Perc_AlCa, Perc_AlMgFe, P(i), Tsolid(i));
        gardensityzero = garnetdensity(Mgnumgar, Perc_AlCa, Perc_AlMgFe, 1e-4, 1);
        gardensitysol = garnetdensity(Mgnumgar, Perc_AlCa, Perc_AlMgFe, 1e-4, Tsolid(i));

        Mgnumalpha = MG(i,2);
        alphadensity = olivinedensity(Mgnumalpha, P(i), Tsolid(i));
        aldensityzero = olivinedensity(Mgnumalpha, 1e-4, 1);
        aldensitysol = olivinedensity(Mgnumalpha, 1e-4, Tsolid(i));

        density = 0.30*cpxdensity + 0.20*gardensity + 0.50*alphadensity;
        densityzero = 0.30*cpxdensityzero + 0.20*gardensityzero + 0.50*aldensityzero;
        densitysol = 0.30*cpxdensitysol + 0.20*gardensitysol + 0.50*aldensitysol;

        D(i) = density;
        Do(i) = densityzero;
        Dsol(i) = densitysol;
    end
end


[initDsolinv, initIndex] = sort(Dsol);	% sorts density from lightest = 1 to densest = endi
for i = 1:endi	% resorts into heaviest at i = 1 and lightest at i = endi, surface
    Dsolinv(i) = initDsolinv(endi+1-i);
    Index(i) = initIndex(endi+1-i);	% "Index" contains the layer number that the density value
end									% was in before being sorted


rinv(1) = R - DM;
Solidus(1) = solidus(RtoP(rinv(1)));
%(-1.1601e-007)*(rinv(1)/1000)^3 + 0.0014*(rinv(1)/1000)^2 + -6.3821*(rinv(1)/1000) + 1.4439e+004 - 200;% melting solidus 200 degrees higher than bulk mantle

for j = 2:endi	% calculate new radius for each layer based on crystalvolume before sorting
    rinv(j) = ((3*(0.001*Mantlevolume))/(4*pi) + ((r(j-1))^3))^(1/3);
    Solidus(j) = solidus(RtoP(rinv(j)));
    % (-1.1601e-007)*(rinv(j)/1000)^3 + 0.0014*(rinv(j)/1000)^2 + -6.3821*(rinv(j)/1000) + 1.4439e+004 - 200;% melting solidus 200 degrees higher than bulk mantle
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for j = 2:endi
%        figure(38); hold on;
%        plot(r(Index(j))/1000,rinv(j)/1000);
% end
% figure(38)
% plot(r/1000, r/1000, 'r:');
% xlabel('Radius before overturn');
% ylabel('Radius after overturn');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hnewatm = 0;
Cnewatm = 0;
for j = 1:endi	% correct temperature for adiabatic rise or fall, and make overturned solid vector
    deltaz(j) = (rinv(j) - r(Index(j)));    % [m]
    dTdZ = (g*alpha*Tsolid(Index(j)))/Cp;   % [C/m]
    Tinv(j) = Tsolid(Index(j));             % overturned temperature without adiabatic correction
    CorrectedTinv(j) = Tsolid(Index(j)) - deltaz(j)*dTdZ; % temp minus adiabatic gradient
    if (CorrectedTinv(j)- Solidus(j))>0
        Tempabovesolidus(j) = CorrectedTinv(j)- Solidus(j);
    else
        Tempabovesolidus(j) = 0;
    end
    test = rinv(j);
    if test < 3200
        meltpercentbottom(j) = Tempabovesolidus(j)*dFdT;
        meltpercenttop(j) = 0;
        Hnewatm = Hnewatm + (0.001*Mantlevolume)*solid(j,10)*0.99;
        Cnewatm = Cnewatm + (0.001*Mantlevolume)*solid(j,10)*0.99;
    else
        meltpercenttop(j) = Tempabovesolidus(j)*dFdT;
        meltpercentbottom(j) = 0;
        Hnewatm = Hnewatm + (0.001*Mantlevolume)*solid(j,10)*0.99;
        Cnewatm = Cnewatm + (0.001*Mantlevolume)*solid(j,10)*0.99;
    end
    meltvolumetop(j) = ((0.001*Mantlevolume).*(meltpercenttop(j))/100);         % in m^3
    meltvolumebottom(j) = ((0.001*Mantlevolume).*(meltpercentbottom(j))/100);   % in m^3
    meltvolume(j) = meltvolumebottom(j) + meltvolumetop(j);                     % in m^3

    Dinv(j) = D(Index(j));
    Doinv(j) = Do(Index(j));
    solidinv(j,:) = solid(Index(j),:);
    SmNdinv(j) = (solid(Index(j),8)*0.15)./(solid(Index(j),9)*0.24); % 147Sm is 15% of Sm; 144Nd is 24% of Nd; 143Nd is 12%
    MGinv(j,:) = MG(Index(j),:);
    LuHf(j)=(solidinv(j,8))/(solidinv(j,9));
    Lu6Hf7inv(j)=0.14*(solidinv(j,8))/(solidinv(j,9)); %(0.0259*solidinv(j,8)*10^4/0.0243)/(0.00162*solidinv(j,9)*10^4/0.1040);
    solidOHinv(j) = solidOH(Index(j));
    intliqOHinv(j) = intliqOH(Index(j));
end

totalmelt = sum(meltvolume(:));  %m^3
%depth = R/1000 - (((R/1000)^3) - (3*totalmelt)/(4*pi))^(1/3);   % [km]
depth = totalmelt/(4*pi*(R)^2)/1000;   % [km]

PEfin=0;PEdelta=0;PEdeltaf=0;PEinit=0;
for i = 1:(endi-1)  % calculate potential energy of final stratigraphy
    PEdeltaf = [Doinv(i)*((4/3)*pi*((rinv(i+1))^3 - (rinv(i))^3))]*g*([(rinv(i+1))+(rinv(i))]/2);
    PEfin = PEfin +PEdeltaf;
    PEdelta = [Do(i)*((4/3)*pi*((r(i+1))^3 - (r(i))^3))]*g*([(r(i+1))+(r(i))]/2);
    PEinit = PEinit +PEdelta;
end
Tempinc = ((PEinit - PEfin)/Cp)/(Planetmass - Coremass); %calculate degrees overturn could raise entire magma ocean
% disp(['Ratio of final to initial potential energy of mantle: ', num2str(PEfin/PEinit)]);
% disp(['Entire mantle temp raised by: ', num2str(Tempinc)]);


% this part to output text files for input to CitCom
% prof = zeros(i+1,2);
% fidr = fopen('radius.txt','w');
% fidp = fopen('density.txt','w');
% fidt = fopen('temp.txt','w');
% countr = fprintf(fidr,'%4.0f\r',r(:));
% countp = fprintf(fidp,'%4.0f\r',Do(:));
% countt = fprintf(fidt,'%4.0f\r',T(:));
% cl = fclose('all');
