a = [1 2 3];
b = [1 2 3];
c = [1 2 3; 4 5 6; 7 8 9];
figure(1);
pcolor(a,b,c);

figure(2);
pcolor(a,b,c); shading interp; colorbar;