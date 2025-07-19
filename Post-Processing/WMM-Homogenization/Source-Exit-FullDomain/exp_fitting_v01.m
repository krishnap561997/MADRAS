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

function [lambda] = exp_fitting_v01(x,y,xlab,ylab,fpath,unique_dp,Vs,ACH,H,fontsize)

for i=1:numel(y(:,1))
    CM = jet(i);
end
lambda = zeros(size(y(:,1)));

figure
%frame_h = get(handle(gcf),'JavaFrame');
%set(frame_h,'Maximized',1);
for i=1:numel(y(:,1))
    g = fittype('1-exp(-a*x)');
%     f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
    f0 = fit(x',y(i,:)',g,'StartPoint',[0.02]);
    MyCoeffs = coeffvalues(f0);

    lambda(i,1) = MyCoeffs(1);
%     xx = linspace(1,8,50);
    
    hold on
    plot(x',y(i,:)','o',x,f0(x),'-','LineWidth',2,'MarkerSize',5,'color',CM(i,:))
    text = strcat(['$',xlab,'$']);
    xlabel(text,'Interpreter','latex')
    text = strcat(['$',ylab,'$']);
    ylabel(text,'Interpreter','latex')
%     text = strcat(['Setup ',num2str(geom_setup)]);
%     title(text)
    set(gca,'FontSize',fontsize)
end
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 16, 9], 'PaperUnits', 'Inches', 'PaperSize', [16 9])

name = strcat(['Fit.pdf']);
print('-dpdf','-r600',fullfile(fpath, name))
name = strcat(['Fit.fig']);
saveas(gcf,fullfile(fpath, name))

figure
%frame_h = get(handle(gcf),'JavaFrame');
%set(frame_h,'Maximized',1);
p1=plot(unique_dp*1e6,lambda,'x','MarkerSize',20,'LineWidth',10);
hold on
p2=plot(unique_dp*1e6,ACH/3600+Vs./H,'o','MarkerSize',20,'LineWidth',10);
title(['Entire room with ACH = ',num2str(ACH)])
xlabel('diameter (\mum)')
% ylabel('\lambda_{c,sim}')
set(gca,'FontSize',fontsize)
h = [p1(1);p2];
text1 = strcat(['\lambda_{c,sim} ; ACH = ',' ',num2str(ACH)]);
text2 = strcat(['\lambda_{c,th}  ; ACH = ',' ',num2str(ACH)]);
legend(h,'\lambda_{c,sim}','\lambda_{c,th}','Location','northwest');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 16, 9], 'PaperUnits', 'Inches', 'PaperSize', [16 9])

name = strcat(['Lambda.pdf']);
print('-dpdf','-r600',fullfile(fpath, name))
name = strcat(['Lambda.fig']);
saveas(gcf,fullfile(fpath, name))

