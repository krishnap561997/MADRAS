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

function [id,dp,tm,tm2,surf] = record_times(id,dp,tm,nfile)
tm2  = 0*tm;
surf = 0*id;

for n=2:nfile
    tic

    infile = strcat(sprintf('parout%05d.h5\n', n));
    disp('-----------------------------------')
    disp(['Reading ',infile])

    outfile    = strcat(sprintf('parout%05d.h5', n));
    data_time  = h5read(outfile,'/g1/time');
    data_id    = h5read(outfile,'/g1/tag');
    data_dp    = h5read(outfile,'/g1/dp');
    data_loc   = h5read(outfile,'/g1/loc');
    data_cycle = h5read(outfile,'/g1/cyc');
    data_wall  = h5read(outfile,'/g1/wll');
    data_vel   = h5read(outfile,'/g1/vel'); % Particle velocity from the solution array (only useful for dp>10mu, it is 0 otherwise)
    data_uel   = h5read(outfile,'/g1/uel'); % Fluid velocity at marker
    data_uelp  = h5read(outfile,'/g1/uelp'); % Langevin perturbation at marker
    data_tot   = h5read(outfile,'/g1/utot'); % Fluid + Langevin at marker
    data_lan   = h5read(outfile,'/g1/lan'); % Some langevin related data

    toc
    fprintf('\n')
% --------------------------------------------------------------------
%% Record the exit time of particles for each output
    tic
    [tf,idx] = ismember(id,data_id);
    tm2(tf) = data_time;
    surf(tf) = data_wall(idx(idx~=0));
    toc
    fprintf('\n')
% --------------------------------------------------------------------
end

