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

clc, clear all, close all;

%% Add Path and Load Input File. Contains par, parout file locations (NOT NEEDED FOR THIS),
clear all; close all; clc;
%location to store data, id_dp_tm location, file size, tlim, room
%dimensions, sphere diameter, distance between source and sink
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f','Delimiter','\n');

fpath   = char(strip(strip(C{1},'left',''''),'right',''''))     
%File dump location. CHANGE IT!
sph     = char(strip(strip(C{2},'left',''''),'right',''''))     %id_dp_tm Location
fclose(fileID);

addpath(genpath(sph))

nfile           = C{3}
time_limit      = C{4}
nxi             = C{5}      %Source Data
nyi             = C{6}
nzi             = C{7}
rlenx           = C{8}
rleny           = C{9}
rlenz           = C{10}
sphere_diameter = C{11}
nxs           = C{12}
nys           = C{13}
nzs           = C{14}
 

%% Load Source Data from id_dp_tm and Making it unique
tic
disp('-----------------------------------')
disp('Load Source data')
[id_source,dp_source,tm_source] = load_source_data(nxi,nyi,nzi,nfile,time_limit);
toc
fprintf('\n')

tic
disp('-----------------------------------')
disp('Make data unique')
[src_id_unique,src_dp_unique,src_tm_unique] = make_unique(id_source,dp_source,tm_source);
toc
fprintf('\n')

%% Determining the locations of all Sinks at the chosen distance from Source
tic
disp('-----------------------------------')
disp('Locating Neighbours')
%Sinks loc for one particular row
toc
fprintf('\n')

%% Cumulative variables
%% Iterating over every Source Sink combination for given source
    
    %%Load Sink Data
    tic
    disp('-----------------------------------')
    disp('Load Sink data')
    nx_snk = nxs;
    ny_snk = nys;
    nz_snk = nzs;
    [snk_id,snk_dp,snk_tm] = load_sink_data(nx_snk,ny_snk,nz_snk);
    toc
    fprintf('\n')
    
    %%Identify particles within Source that enter Sink
    tic
    disp('-----------------------------------')
    disp('Comparing time taken for particles to reach from Source to Sink')
    [src2snkDt] = src2snk_v02(src_id_unique,src_tm_unique,snk_id,snk_tm);
    %[snk_id_unique, snk_dp_unique, src2snkDt_unique] = make_unique(snk_id, snk_dp, src2snkDt);
    toc
    
    %%Store and Plot Source to Sink dt
    tic
    disp('-----------------------------------')
    disp('Comparing time taken for particles to reach from Source to Sink')
    [Nexit,Nhat,Ntot] = src2snkstats_v01(src_dp_unique,snk_dp,src2snkDt,time_limit);
    

%% Save file
    tic
    disp('-----------------------------------')
    disp('Saving the workspace')
    text = strcat(['Sink_nxi_',sprintf('%03d',round(nxs)),...
            '_nyi_',sprintf('%03d',round(nys)),...
	    '_nzi_',sprintf('%03d',round(nzs))]);
    source_loc = [nxi, nyi, nzi];
    save(fullfile(fpath, text), 'src2snkDt', 'Nexit',...
        'Ntot', 'Nhat')
    disp('Workspace Saved')
    toc

a+b+c
