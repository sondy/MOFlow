% eer_graph.m

figure(44); %title(['Reference density with depth pre- and post-overturn for model:  ', name]); 
    hold on; 
    plot(Dsolinv, 1:j, 'k')%, 'LineWidth', 4); %        plot(Dsol(1:marker2), r(1:marker2)/1000,'y',...     
    
    d1 = diff(Dsolinv);
    
    d2 = diff(d1);
    
    plot(2:j, d1, 'r')
    
    plot(3:j, d2, 'b')
    
    xlabel('density at 1 atm and solidus temperature [kg m^{-3}]'); 
    ylabel('radius, km'); 
%    xlim([2600 4200])
%    ylim([2500 6500])
    
    %legend('Pre-overturn density', 'Post-overturn density', 'Location', 'East')
    
    justDensityWithDepth = strcat('plots/justDensityWithDepth',...
        DM_string, '.eps');
    print('-depsc', justDensityWithDepth)