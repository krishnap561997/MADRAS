%% *** MAtlab turbulent Dispersion Room-scale Analysis Software (MADRAS)*** 
% 
% This software is for calculating turbulent dispersal statistics, such as
% concentration and flux, from Lagrangian particle tracks generated from
% DNS, LES, and RANS.
% 
% The software has also been specialized for predicting pathogen concentrat
% -ion in indoor spaces to compute quantities such as cumulative exposure
% time (CET) and Safe Occupancy Limit
% 
% The use/distribution of this code is allowed ONLY after approval by Prof.
% S. Balachandar *(bala1s@ufl.edu)*, the project PI. This code was written
% and conceptualized at the University of Florida at the Computational
% Multiphysics Group and the Center for Compressible Multiphase Turbulence
% by Prof. S. Balachandar, Prof. Nadim Zgheib, Dr. Jorge Salinas, and K.A.
% Krishnaprasad.
% *************************************************************************

function [] = figure_v01(CM,unique_dp,xval,yval,xlab,ylab,titlelab,fpath,fontsize)

resolution = 300; % dpi

figure
for ii=1:numel(unique_dp)
    plot(xval,yval(ii,:),'-','LineWidth',5,'color',CM(ii,:))
    hold on

    Legend{ii} = strcat(['d_p =',' ',num2str(unique_dp(ii)*1e6),' \mum']);   
    legend(Legend,'Location', 'eastoutside');
    text = strcat(['$',xlab,'$']);
    xlabel(text,'Interpreter','latex')
    text = strcat(['$',ylab,'$']);
    ylabel(text,'Interpreter','latex')
%     text = strcat(['Setup ',num2str(geom_setup)]);
%     title(text)
    set(gca,'FontSize',fontsize)
end
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 16, 9], 'PaperUnits', 'Inches', 'PaperSize', [16 9])

name = strcat([titlelab,'.pdf']);
print('-dpdf','-r600',fullfile(fpath, name))
% saveas(gcf,name)
% exportgraphics(gcf,namepng,'Resolution',resolution)
name = strcat([titlelab,'.fig']);
saveas(gcf,fullfile(fpath, name))

name = strcat([titlelab,'.png']);
print('-dpng','-r600',fullfile(fpath, name))
