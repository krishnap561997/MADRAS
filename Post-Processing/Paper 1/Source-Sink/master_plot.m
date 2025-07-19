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

clear all; close all; clc;
path =  '/blue/bala/krishnap.kalivel/COVID/room_29/Source2Sink/Data_2m'
addpath(genpath(path))

fpath = '/blue/bala/krishnap.kalivel/COVID/Sou2Snk_Spheres/dist_2m/room_29'

fontsize = 18;  nxt = 25;   nyt = 8;    nzt = 25; time_limit = 700;
comb_n = 87918;
d = 0.4;
V_sink = pi*(d^2)/4;
unique_dp = [0.2 1 5 10 15 20 25 30 40 50];
%Nexit_c =zeros(10, 701, comb_n);
%Nhat_c =zeros(10, 701,comb_n);
%Ntot_c = zeros(10,1, comb_n);
counter = 1;
%sinks_loc_c = zeros(comb_n,3);
%source_loc_c =zeros(comb_n,3);
for nxi = 1:nxt
    for nyi = 1:nyt
        for nzi = 1:nzt
            [nxi nyi nzi]
            name = strcat(['Source_nxi_',sprintf('%03d',round(nxi)),...
                                  '_nyi_',sprintf('%03d',round(nyi)),...
                                  '_nzi_',sprintf('%03d',round(nzi)),'.mat']);
            load(name)
            for j=1:numel(Nexit_cum)
                    Nexit_c(:, :, counter) = cell2mat(Nexit_cum(j));
                    Nhat_c(:, :, counter) = cell2mat(Nhat_cum(j));
                    Ntot_c(:,:,counter) = cell2mat(Ntot_cum(j));
                    sinks_loc_c(counter,:) =sinks_loc(j,:);
                    source_loc_c(counter,:) = source_loc;
                    counter = counter + 1;
            end
        end
    end
end


%Plotting the mean of Nhat and one standard deviation above and below the mean
%Nhat_mean = mean(Nhat_c,3);
%Nhat_std = zeros(size(Nhat_c,1),size(Nhat_c,2));

%for i=1:size(Nhat_c,1)
%	for j=1:size(Nhat_c,2)
%		Nhat_std(i,j) = std(Nhat_c(i,j,:));
%	end
%end
%
%Nhat_ms1 = Nhat_mean + Nhat_std;
%Nhat_ms2 = Nhat_mean - Nhat_std;
%
%for i=1:numel(unique_dp)
%    figure(i);
%    plot(0:time_limit,Nhat_mean(i,:),'-k','LineWidth',2);
%    hold on;
%    plot(0:time_limit, Nhat_ms1(i,:),'-r','LineWidth',2);
%    plot(0:time_limit, Nhat_ms2(i,:),'-b','LineWidth',2);
%    ylabel('$\widehat{N}$','Interpreter','latex')
%    xlabel('t (s)','Interpreter','latex')
%    set(gca,'FontSize',18)
%    title({strcat('dp: ',num2str(unique_dp(i)),'$\mu$m')}, 'FontSize', 16, 'Interpreter', 'latex')
%    legend({'$\overline{\widehat{N}}$', '$\overline{\widehat{N}}$ + $\sigma$',...
%        '$\overline{\widehat{N}}$ - $\sigma$'},'Interpreter', 'latex', 'Location', 'southeast','FontSize',16);
%    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 16, 9], 'PaperUnits', 'Inches', 'PaperSize', [16 9])
%    name = strcat('Nhat_mean_std_',num2str(unique_dp(i)),'micron.pdf');
%    print('-dpdf','-r300',fullfile(fpath, name))
%    name = strcat('Nhat_mean_std_',num2str(unique_dp(i)),'micron.fig');
%    saveas(gcf,fullfile(fpath, name))
%    name = strcat('Nhat_mean_std_',num2str(unique_dp(i)),'micron.png');
%    print('-dpng','-r300',fullfile(fpath, name))	
%end



%close all;

Nexit_mean = mean(Nexit_c,3);
Nexit_std = zeros(size(Nexit_c,1),size(Nexit_c,2));
for i=1:size(Nexit_c,1)
	for j=1:size(Nexit_c,2)
		Nexit_std(i,j) = std(Nexit_c(i,j,:));
	end
end

Nexit_ms1 = Nexit_mean + Nexit_std;
Nexit_ms2 = Nexit_mean - Nexit_std;

close all;
for i=1:numel(unique_dp)
    figure(i);
    bar(0:time_limit,Nexit_mean(i,:), 'facecolor','flat');
    hold on;
    plot(0:time_limit,Nexit_ms1(i,:), '--r','LineWidth', 2);
    plot(0:time_limit,Nexit_ms2(i,:), '--b','LineWidth', 2);
    ylabel('$\overline{N}$','Interpreter','latex')
    xlabel('t (s)','Interpreter','latex')
    set(gca,'FontSize',18)
    legend({'$\overline{N}$', '$\overline{N}$ + $\sigma$',...
        '$\overline{N}$ - $\sigma$'},'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 16);
    title({strcat('dp: ',num2str(unique_dp(i)),'$\mu$m')}, 'FontSize', 16, 'Interpreter', 'latex')
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 16, 9], 'PaperUnits', 'Inches', 'PaperSize', [16 9])
    name = strcat('Nexit_mean_',num2str(unique_dp(i)),'micron.pdf');
    print('-dpdf','-r300',fullfile(fpath, name))
    name = strcat('Nexit_mean_',num2str(unique_dp(i)),'micron.fig');
    saveas(gcf,fullfile(fpath, name))
    name = strcat('Nexit_mean_',num2str(unique_dp(i)),'micron.png');
    print('-dpng','-r300',fullfile(fpath, name))
end

