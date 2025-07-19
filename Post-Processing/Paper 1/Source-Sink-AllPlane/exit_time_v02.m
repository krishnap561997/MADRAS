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

function [Nexit,Nhat,Ntot] = src2snkstats_v01(src_dp_unique,snk_dp,snk_tm,src2snkDt,time_limit)

%% Compute the number of particle sizes involved
for i=1:numel(src_dp)
    unique_dp{i,1}   = unique(src_dp{i});
end
Ntot    = zeros(numel(src_dp),numel(unique_dp{1}));
Nexit   = zeros(numel(src_dp),numel(unique_dp{1}),time_limit+1,...
               numel(snk_dp));
Nhat    = zeros(numel(src_dp),numel(unique_dp{1}),time_limit+1,...
               numel(snk_dp));
