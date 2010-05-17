% cool2clementwhole.m, FOR WHOLE-MANTLE MAGMA OCEANS

disp('In cool2clementwhole')

% program to calculate temperature evolution of the Martian mantle, with heat production
% bottom of magma ocean is defined by an end point node - core temperature allow to evolve
% through time
% for whole-mantle MOs change lines marked with XXX

Rcore = (R - DM)/1000;		% km, bottom of MO in this case not CMB
Rmantle = R/1000;           % km
tempsurf = Tinv(j) + 273;   % surface temp from end of sortandinvert
rhomantle = solidrho;       % kg/m^3

nx = 990;                   % number of nodes; 990 as set in MOFlow XXXXXXXXXXXXXXX
dx =(Rmantle-Rcore)/nx;     % km per node
dt = (dx*1e3)^2/(2*kappa);  % sec
nt = round(tfinal/dt);

time2c = zeros(nt+1,1);
r2c = zeros(nx+2,1);
temp = zeros(nx+2,nt+1);
xchond = zeros(nx+2,1);
newtemp = zeros(nx+2,1);
quartertemp = zeros(nx+2,1);
halftemp = zeros(nx+2,1);
threeqtemp = zeros(nx+2,1);
coreflux = zeros(nt+1,1);
surfaceflux = zeros(nt+1,1);
shallowsolidus = zeros(nx+2,1);
CMBsolidus = zeros(nx+2,1);

% initialization 

r2c(1) = Rcore-dx;
temp(1,1) = tempcore;
for ix = 2:990+1   % creating temperature profile from sortandinvert.m
    r2c(ix,1) = Rcore+(ix-2)*dx;
    temp(ix,1) = Tinv(ix-1) + 273;         
end
temp(992) = temp(991); 
r2c(992) = R/1000 + dx;

counter = 1;

for it = 0:dt:tfinal-dt
   
   time2c(counter) = (it/(3.14*1e13));
%    Heatu238=hu238*cu238*exp((tage-it)*log(2)/thu238); Heatu235=hu235*cu235*exp((tage-it)*log(2)/thu235);
%    Heatk40=hk40*ck40*exp((tage-it)*log(2)/thk40); Heatth232=hth232*cth232*exp((tage-it)*log(2)/thth232);   
%    Hstep=Heatu238+Heatu235+Heatk40+Heatth232;          % heat for this timestep in W/kg-mantle
   
   newtemp(1) = tempcore;     % these two are boundary conditions for the conductive heat profile
   newtemp(nx+2) = tempsurf; 
   
   for ix = 2:nx+1 % calculate the temperature at each node from core to surface
      term1 = (1/((r2c(ix)*1e3)*(dx*10^3)))*(temp(ix+1,counter)-temp(ix-1,counter));
      term2 = (1/(dx*10^3)^2)*(temp(ix+1,counter)-2*temp(ix,counter)+temp(ix-1,counter));
%       term3=Hstep/Cp*xchond(ix);
      newtemp(ix) = temp(ix,counter)+(kappa*dt)*(term1+term2); %+(dt*term3);
      
      if it == round(nt/4)*dt
          quartertemp(ix) = newtemp(ix);   % for graphing temp profile quarterway
      end
      if it == round(nt/2)*dt
          halftemp(ix) = newtemp(ix);   % for graphing temp profile halfway
      end
      if it == round(3*nt/4)*dt
          threeqtemp(ix) = newtemp(ix);   % for graphing temp profile threequarterway
      end     
      if round(newtemp(ix) - 1273) < 30  % finds height at which temp becomes <=1000C and records for plotting
            if round(newtemp(ix) - 1273) > 0
                crustlevel = r2c(ix);
            end
      end
   end
   
   coreflux(counter) = (kappa*(rhomantle*Cp))*(newtemp(1)-newtemp(2))/(dx*10^3);
   tempcorenew = tempcore-(dt*coreflux(counter)*3/(rhocore*cpcore*(Rcore*1e3)));    % XXX
   tempcore = tempcorenew;
   surfaceflux(counter) = endemiss*sigma*((newtemp(nx+2))^4 - 250^4);
   tempsurfnew = tempsurf-(dt*surfaceflux(counter)*3/(rhomantle*Cp*(Rmantle*1e3)));
   tempsurf = tempsurfnew;
   temp(:,counter+1)=newtemp(:,1);
   
counter = counter+1;   
end
tempfinal = temp(:,counter);

