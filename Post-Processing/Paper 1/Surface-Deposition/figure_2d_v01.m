function [] = figure_2d_v01(soln,xmesh,ymesh,diam,xlab,ylab,titlelab,fpath,fontsize)

h = fspecial('gaussian');
soln_filtered = filter2(h, soln);
for ii=1:0
    soln_filtered = filter2(h, soln_filtered);
end
contourf(xmesh,ymesh,soln_filtered,'edgecolor','none');

colormap(hot(256));
colorbar;
axis on;
axis equal
text = strcat(['$',xlab,'$']);
xlabel(text,'Interpreter','latex')
text = strcat(['$',ylab,'$']);
ylabel(text,'Interpreter','latex')
text = strcat(['$',titlelab,'\,$ ; ','$d_p = $',num2str(diam*1e6),...
    '$\mu m$']);
title(text,'Interpreter','latex')
set(gca,'FontSize',fontsize)

namepng = strcat([titlelab,'_','d_p=',num2str(diam*1e6),'micron.png']);
saveas(gcf, fullfile(fpath, namepng));

namefig = strcat([titlelab,'_','d_p=',num2str(diam*1e6),'micron.fig']);
saveas(gcf, fullfile(fpath, namefig));
            
