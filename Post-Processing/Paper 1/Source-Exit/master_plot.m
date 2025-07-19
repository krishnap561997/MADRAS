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

clc, clear all, close all
path = '/blue/bala/krishnap.kalivel/COVID/room_29/Source2Exit_tlim700/Sou2ExtStats'
addpath(genpath(path))

fpath = '/blue/bala/krishnap.kalivel/COVID/room_29/Source2Exit_tlim700/Figures'
fontsize = 25;  nxt = 25;   nyt = 8;    nzt = 25; time_limit = 700;
Nexit_c   = zeros(10,time_limit+1,nxt*nyt*nzt); %unique_dp,number of spheres,time_limit
Nwall_c   = zeros(10,time_limit+1,nxt*nyt*nzt); %unique_dp,number of spheres,time_limit
Noutlet_c = zeros(10,time_limit+1,nxt*nyt*nzt); %unique_dp,number of spheres,time_limit
Nfloor_c  = zeros(10,time_limit+1,nxt*nyt*nzt); %unique_dp,number of spheres,time_limit
counter = 0;
for nxi = 1:nxt
    for nyi = 1:nyt
        for nzi = 1:nzt
            counter = counter + 1
            [nxi nyi nzi]
            name = strcat(['Sou2Ext_nxi_',sprintf('%03d',round(nxi)),...
                                  '_nyi_',sprintf('%03d',round(nyi)),...
                                  '_nzi_',sprintf('%03d',round(nzi)),'.mat']);
            load(name)
            for i=1:10
                Nexit_c(i,:,counter) = Nexit(i,:)/Ntot(i);
                Nwall_c(i,:,counter) = Nwall(i,:)/Ntot(i);
                Nfloor_c(i,:,counter) = Nfloor(i,:)/Ntot(i);
                Noutlet_c(i,:,counter) = Noutlet(i,:)/Ntot(i);
            end
        end
    end
end
                
%            for i=1:10
%                tic
%                figure(i)
%                subplot(2,2,1)
%                toto(1:time_limit+1,:) = Nexit_c(i,:,:);
%                plot(0:time_limit,toto,'.')
%                xlabel('$t\,(sec)$','Interpreter','latex')
%                ylabel('$N_{exit}$','Interpreter','latex')
%                set(gca,'FontSize',fontsize)
%                subplot(2,2,2)
%                toto(1:time_limit+1,:) = Nwall_c(i,:,:);
%                plot(0:time_limit,toto,'.')
%                xlabel('$t\,(sec)$','Interpreter','latex')
%                ylabel('$N_{wall}$','Interpreter','latex')
%                set(gca,'FontSize',fontsize)
%                subplot(2,2,3)
%                toto(1:time_limit+1,:) = Nfloor_c(i,:,:);
%                plot(0:time_limit,toto,'.')
%                xlabel('$t\,(sec)$','Interpreter','latex')
%                ylabel('$N_{floor}$','Interpreter','latex')
%                set(gca,'FontSize',fontsize)
%                subplot(2,2,4)
 %               toto(1:time_limit+1,:) = Noutlet_c(i,:,:);
%                plot(0:time_limit,toto,'.')
%                xlabel('$t\,(sec)$','Interpreter','latex')
%                ylabel('$N_{outlet}$','Interpreter','latex')
%                set(gca,'FontSize',fontsize)
%                toc
%                tic
%                set(gcf, 'Units', 'Inches', 'Position', [0, 0, 16, 9], 'PaperUnits', 'Inches', 'PaperSize', [16 9])
%                name = strcat(['dp_',num2str(unique_dp(i)*1e6),'micron.png']);
%                print('-dpng','-r300',fullfile(fpath, name))
%                name = strcat(['dp_',num2str(unique_dp(i)*1e6),'micron.fig']);
%                saveas(gcf,fullfile(fpath, name))
 %               toc
%            end
text = 'Data_tlim700';
save(fullfile(fpath, text), 'Noutlet_c','Nexit_c','Nfloor_c','Nwall_c');
a+b+c
