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

function [Nexit,Nwall,Nfloor,Noutlet,Ntot,unique_dp] ...
    = exit_time(dp,time,time2,nfile,time_limit,exit_surf)

%% Compute the number of particle sizes involved
unique_dp = unique(dp);
    
Ntot    = zeros(numel(unique_dp),1);
Nwall   = zeros(numel(unique_dp),time_limit+1);
Nfloor  = zeros(numel(unique_dp),time_limit+1);
Noutlet = zeros(numel(unique_dp),time_limit+1);
Nexit   = zeros(numel(unique_dp),time_limit+1);

for j=1:numel(unique_dp) % Number of particle sizes
    idx         = dp == unique_dp(j); 
    Ntot(j)   = sum(idx);
    for k = 1:time_limit+1 % Number of time instances
        idx = (dp == unique_dp(j) ...
            & (time2 > time) ...
            & (time2 - time < k-1)...
            & exit_surf>1 &  exit_surf<7);
        Nwall(j,k)  = sum(idx);
        
        idx = (dp == unique_dp(j) ...
            & (time2 > time) ...
            & (time2 - time < k-1)...
            & exit_surf==1);
        Nfloor(j,k)  = sum(idx);
        
        idx = (dp == unique_dp(j) ...
            & (time2 > time) ...
            & (time2 - time < k-1)...
            & exit_surf==7);
        Noutlet(j,k)  = sum(idx);

        Nexit(j,k)  = Noutlet(j,k) + Nfloor(j,k) + Nwall(j,k);
    end
end

