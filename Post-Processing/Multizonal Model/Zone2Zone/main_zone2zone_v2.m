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
%% Initial particles
% This routine computes the residence time of particles for entire room

%% Add path
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %s %f %f %f %f %f %f %f %f %f','Delimiter','\n');

fpath   = char(strip(strip(C{1},'left',''''),'right',''''))
spath  = char(strip(strip(C{2},'left',''''),'right',''''))
fclose(fileID);

nfile    = C{3}
rlenx    = C{4}
rleny    = C{5}
rlenz    = C{6}
outlet_x = C{7}
outlet_z = C{8}
start_t  = C{9}
Z_src    = C{10}
dw       = C{11}

outlet_c = [outlet_x, outlet_z];
%% Load Source Zone Data
tic
disp('-----------------------------------')
disp(strcat(['Loading Source Data at t = ',num2str(start_t),' s from zone ',num2str(Z_src)]));
[id_src0, dp_src0, loc_src0] = load_srczone_data_v2(fpath, start_t, Z_src, rlenx, rleny, rlenz, outlet_c, dw);
toc
fprintf('\n')

%% Source to Sink Statistics
time = 0:(nfile - start_t - 1);
unique_dp = unique(dp_src0);

for i=1:numel(unique_dp)
    id_src{i} = id_src0(dp_src0 == unique_dp(i));
    Ntot(i) = numel(id_src{i});
end

Nhat = zeros(numel(unique_dp),numel(time),4);
Nexit = zeros(numel(unique_dp),numel(time),4);


tic
disp('-----------------------------------')
disp('Reading zonal data at sink')
for i=2:(nfile-start_t)
    count = get_zonal_stats(fpath, start_t, i, rlenx, rlenz, outlet_c, id_src, unique_dp);
    Nexit(:,i,:) = reshape(count,[numel(unique_dp),1,4]);
    Nhat(:,i,:) = Nhat(:,i-1,:) + reshape(count,[numel(unique_dp),1,4]);
end

toc
%% Saving Data
tic
disp('-----------------------------------')
disp('Saving Data')
sname = strcat(['IZData_Zone_',num2str(Z_src),'_t_',sprintf('%04d',start_t),'s.mat']);
save(fullfile(spath,sname),'Nhat','Nexit','Ntot','Z_src','start_t','unique_dp');
