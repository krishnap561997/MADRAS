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
%Average Concentration for ACH 10
%lambda2 = [1.2197    1.2185    1.1803    0.9825    0.7429    0.5499    0.4081    0.3086    0.1870    0.1205];

%Nexit_c =zeros(10, 701, comb_n);
%Nhat_c =zeros(10, 701,comb_n);
%Ntot_c = zeros(10,1, comb_n);
%sinks_loc_c = zeros(comb_n,3);
%source_loc_c =zeros(comb_n,3);
for it = 1:5
    path =  '/blue/bala/krishnap.kalivel/COVID/Sou2Snk_Spheres/dist_2m/room_27/Data';
    
    addpath(genpath(path))
    
    fpath = '/blue/bala/krishnap.kalivel/COVID/Sou2Snk_Spheres/dist_2m/room_27/Figures';
        
    s = (it-1)*5 +1;
    e = it*5;
    
    counter = 1;
    unique_dp = [0.2 1 5 10 15 20 25 30 40 50];
    N_dp = numel(unique_dp);
    fontsize = 20;  nxt = 25;   nyt = 8;    nzt = 25; time_limit = 1900;
    d = 0.4;
    V_sink = (4/3)*pi*(d^3);
    nxt1 = 5;

    for nxi = s:e
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
                        Ntot_c = cell2mat(Ntot_cum(j)).*ones(size(Nhat_c,1),size(Nhat_c,2));
                        Nexit_c(:,:,counter) = Nexit_c(:,:,counter)./Ntot_c;
                        Nhat_c(:,:,counter) = Nhat_c(:,:,counter)./Ntot_c;
                        sinks_loc_c(counter,:) =sinks_loc(j,:);
                        source_loc_c(counter,:) = source_loc;
                        counter = counter + 1;
                end
            end
        end
    end
    comb_n = size(Nexit_c,3);
    Nhat_end = reshape(Nhat_c(:,end,:), [N_dp, comb_n]);
    Nhat_sum = sum(Nhat_c, 3);
    Nexit_sum = sum(Nexit_c,3);
    text = strcat('FigureData_nx_', sprintf('%02d', round(s)),'_',sprintf('%02d', round(e)));
    save(fullfile(fpath,text),'Nhat_end','Nhat_sum', 'Nexit_sum', 'comb_n','sinks_loc_c','source_loc_c');
    clear all; 
end

a+b+c

