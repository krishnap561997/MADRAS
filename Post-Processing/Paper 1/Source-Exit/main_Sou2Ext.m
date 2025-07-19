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
% The routine computes residence time of particles from a finite source

%% Add path
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %s %s %s %f %f %f %f %f','Delimiter','\n');

par     = char(strip(strip(C{1},'left',''''),'right',''''))
parout  = char(strip(strip(C{2},'left',''''),'right',''''))
fpath   = char(strip(strip(C{3},'left',''''),'right',''''))
sph     = char(strip(strip(C{4},'left',''''),'right',''''))
fclose(fileID);

addpath(genpath(par))
addpath(genpath(parout))
addpath(genpath(sph))

nfile           = C{5}
time_limit      = C{6}
nxi             = C{7}
nyi             = C{8}
nzi             = C{9}

tic
disp('-----------------------------------')
disp('Load sphere data')
[id,dp,tm] = load_sphere_data(nxi,nyi,nzi,nfile,time_limit);
toc
fprintf('\n')

tic
disp('-----------------------------------')
disp('Make data unique')
[id_unique,dp_unique,tm_unique] = make_unique(id,dp,tm);
toc
fprintf('\n')

%% Identify particles within the sphere that leave the room
tic
disp('-----------------------------------')
disp('Recording Exit Times')
[id,dp,tm,tm2,surf] = record_times(id_unique,dp_unique,tm_unique,nfile);
toc
fprintf('\n')

%% Verify unaccounted particles are removed
% Unaccounted particles are particles that do not appear in the par files,
% i.e. they have exited the domain, but do ont also appear in the parout
% files.
tic
disp('-----------------------------------')
disp('Removing unaccounted particles')
[id,dp,tm,tm2,surf] = verify_remove(id,dp,tm,tm2,surf,nfile);
toc
fprintf('\n')

%% Generate exit files on walls, floor, ceiling, and outlet
tic
disp('-----------------------------------')
disp('Forming Nexit, Noutlet, Nwall, Nfloor')
[Nexit,Nwall,Nfloor,Noutlet,Ntot,unique_dp] = exit_time(dp,tm,tm2,nfile,time_limit,surf);
toc
fprintf('\n')

%% Save Nexit, Nwall, Nfloor, Noutlet, Ntot
tic
disp('-----------------------------------')
disp('Saving Nexit, Nwall, Nfloor, Noutlet, Ntot')
    text = strcat(['Sou2Ext_nxi_',sprintf('%03d',round(nxi)),...
                          '_nyi_',sprintf('%03d',round(nyi)),...
                          '_nzi_',sprintf('%03d',round(nzi))]);
    save(fullfile(fpath, text),'Nexit','Nwall','Nfloor','Noutlet','Ntot','unique_dp');
    disp('Variables saved')
    toc

disp('-----------------------------------')
disp('Job terminated successfully')
a+b+c % An error message is needed for the job to terminate

