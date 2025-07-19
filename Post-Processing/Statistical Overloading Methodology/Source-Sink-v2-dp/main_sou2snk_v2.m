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
C = textscan(fileID,'%s %s %f %f %f %f %f %f %f %f %f %f %f %f','Delimiter','\n');

fpath   = char(strip(strip(C{1},'left',''''),'right',''''));        %Save location
dpath     = char(strip(strip(C{2},'left',''''),'right',''''));      %id_dp_tm location

nfile           = C{3}      %Number of par files
time_limit_min  = C{4}      %Minimum time limit for Source to Sink Stats
time_gap        = C{5}      %Gap between ejection periods
nxi             = C{6}      %Source x,y,z indices
nyi             = C{7}
nzi             = C{8}
rlenx           = C{9}      %Room size, x, y, z
rleny           = C{10}
rlenz           = C{11}
sphere_diameter = C{12}     %Sink diameter
dist            = C{13}     %Src to Snk distance
tolerance       = C{14}     %Tolerance for distance

%%Load Source Data from id_dp_tm and Making it unique
dx = (sphere_diameter/2):sphere_diameter:(rlenx-sphere_diameter/2);
dy = (sphere_diameter/2):sphere_diameter:(rleny-sphere_diameter/2);
dz = (sphere_diameter/2):sphere_diameter:(rlenz-sphere_diameter/2);
xc_src = dx(nxi); yc_src = dy(nyi); zc_src = dz(nzi);

time_start = 0:(time_gap):(nfile-time_limit_min-1-time_gap);
time_end = time_start + time_gap;
time_limit = max(time_limit_min.*ones(size(time_start)), nfile-time_end-1);

tic
disp('-----------------------------------')
disp('Load Source data')
[id_source,dp_source,tm_source,~,~,~] = load_source_data_v2(nxi,nyi,nzi,dpath);
toc
fprintf('\n')

%%Determining the locations of all Sinks at the chosen distance from Source
tic
disp('-----------------------------------')
disp('Locating Neighbours')
[sinks_loc] = neighbours(nxi, nyi, nzi, rlenx, rleny, rlenz, dist, sphere_diameter, tolerance);
%Sinks loc is a NX3 matrix with x,y,z location of all sinks at a distance
%dist from source
toc
fprintf('\n')

%%Cumulative variables
Nexit_cum       = cell(size(sinks_loc,1),numel(time_start));
Ntot_cum        = zeros(size(sinks_loc,1),numel(time_start));
Nhat_cum        = cell(size(sinks_loc,1),numel(time_start));

%%Iterating over every Source Sink combination for given source
for i=1:size(sinks_loc,1)
    
    %%Load Sink data
    tic
    nx_snk = sinks_loc(i,1);
    ny_snk = sinks_loc(i,2);
    nz_snk = sinks_loc(i,3);
    text = strcat(['Sink: nxi ',num2str(nx_snk),...
        ', nyi ',num2str(ny_snk),...
        ', nzi ',num2str(nz_snk)]);
    disp(text);
    tic
    disp('-----------------------------------')
    disp('Load Sink data')

    [snk_id,snk_dp,snk_tm] = load_sink_data_v2(nx_snk,ny_snk,nz_snk,dpath);
    toc
    fprintf('\n')


    for j=1:numel(time_start)
        tic
        disp('-----------------------------------')
        disp('Remove data from t = 0 to time_gap')
        [id_source0,dp_source0,tm_source0] = remove_data_till_tm(id_source,dp_source,...
            tm_source,time_start(j));
        toc
        fprintf('\n')

        tic
        disp('-----------------------------------')
        disp('Make data unique')
        [src_id_unique0,src_dp_unique0,src_tm_unique0] = make_unique(id_source0,dp_source0,tm_source0);
        toc
        fprintf('\n')

        %%Remove source particles after nfile-time_limit
        tic
        disp('-----------------------------------')
        disp('Removing data from t = time_limit to t = nfile -1')
        [src_id_unique,src_dp_unique,src_tm_unique] = remove_data_after_tm(src_id_unique0,src_dp_unique0,src_tm_unique0,time_end(j));
        toc
        
        %%Identify particles within Source that enter Sink
        tic
        disp('-----------------------------------')
        disp('Comparing time taken for particles to reach from Source to Sink')
        [src2snkDt] = src2snk_v02(src_id_unique,src_tm_unique,snk_id,snk_tm);
        [Nexit_cum{i,j},Nhat_cum{i,j},Ntot_cum(i,j)] = src2snkstats_v02(src_dp_unique,src2snkDt,time_limit(j));
        toc
    end
end

%% Saving workspace
tic
disp('-----------------------------------')
disp('Saving the workspace')
text = strcat(['Source_nxi_',sprintf('%03d',round(nxi)),...
        '_nyi_',sprintf('%03d',round(nyi)),...
    '_nzi_',sprintf('%03d',round(nzi))]);
source_loc = [nxi, nyi, nzi];
save(fullfile(fpath, text), 'Nexit_cum',...
    'Ntot_cum', 'Nhat_cum','sinks_loc','time_start','time_end','time_limit')
disp('Workspace Saved')
toc
