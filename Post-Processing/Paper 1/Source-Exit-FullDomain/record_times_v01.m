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

function [id_1,dp_1,time,exit_surf,exit_loc] = record_times_v01(nfile)

n = 1;
tic
infile = strcat(sprintf('par%05d.h5\n', n));
disp('-----------------------------------')
disp(['Reading ',infile])
time_1  = h5read(infile,'/g1/time');% single
id_1    = h5read(infile,'/g1/tag'); % unit32
dp_1    = h5read(infile,'/g1/dp');  % single
loc_1   = h5read(infile,'/g1/loc'); % single
toc
fprintf('\n')
% --------------------------------------------------------------------
time = 0*single(id_1);
exit_surf = 0*single(id_1);
exit_loc  = 0*single(loc_1);

for n=2:nfile
    tic

    infile = strcat(sprintf('parout%05d.h5\n', n));
    disp('-----------------------------------')
    disp(['Reading ',infile])

    outfile    = strcat(sprintf('parout%05d.h5', n));
    data_time  = h5read(outfile,'/g1/time');
    data_id    = h5read(outfile,'/g1/tag');
%    data_dp    = h5read(outfile,'/g1/dp');
    data_loc   = h5read(outfile,'/g1/loc');
%    data_cycle = h5read(outfile,'/g1/cyc');
    data_wall  = h5read(outfile,'/g1/wll');
%    data_vel   = h5read(outfile,'/g1/vel'); % Particle velocity from the solution array (only useful for dp>10mu, it is 0 otherwise)
%    data_uel   = h5read(outfile,'/g1/uel'); % Fluid velocity at marker
%    data_uelp  = h5read(outfile,'/g1/uelp'); % Langevin perturbation at marker
%    data_tot   = h5read(outfile,'/g1/utot'); % Fluid + Langevin at marker
%    data_lan   = h5read(outfile,'/g1/lan'); % Some langevin related data
    
    toc
    fprintf('\n')
% --------------------------------------------------------------------
%% Record the exit time of particles for each output
    tic
    [tf,idx] = ismember(id_1,data_id);
    exit_surf(tf) = data_wall(idx(idx~=0));  
    exit_loc(tf,:) = data_loc(idx(idx~=0),:);  
    time(tf) = data_time;

    disp([num2str(numel(data_id)),' particles have exited the domain'])
    toc
    fprintf('\n')
% --------------------------------------------------------------------
end

% --------------------------------------------------------------------
%% Use below structure for old parout format in room_24_nolan
% --------------------------------------------------------------------

%for n=2:nfile
%    tic
%    infile = strcat(sprintf('parout%05d.h5\n', n));
%    disp('-----------------------------------')
%    disp(['Reading ',infile])
%    outfile  = strcat(sprintf('parout%05d.h5', n));
%    data_time  = h5read(outfile,'/g1/time');
%    data_id  = h5read(outfile,'/g1/tag');
%    data_dp  = h5read(outfile,'/g1/dp');
%    data_loc = h5read(outfile,'/g1/loc');
%    data_cycle = h5read(outfile,'/g1/cyc');
%    data_wall = h5read(outfile,'/g1/wll');
%    toc
%    fprintf('\n')
%% --------------------------------------------------------------------
%%% Record the exit time of particles for each output
%    tic
%    [tf,idx] = ismember(id_1,data_id);
%    exit_surf(tf) = data_wall(idx(idx~=0));  
%    time(tf) = data_time;
%
%    disp([num2str(numel(data_id)),' particles have exited the domain'])
%    toc
%    fprintf('\n')
%% --------------------------------------------------------------------
%end


