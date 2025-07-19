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

function [Nexit,Nhat,Ntot] = src2snkstats_v02(src_dp_unique, src2snkDt,time_limit)

Ntot = numel(src_dp_unique);
Nexit   = zeros(1,time_limit+1);
Nhat    = zeros(1,time_limit+1);
%% Compute the number of particle sizes involved
for j=1:time_limit+1
    idx = (src2snkDt == j-1);
    Nexit(1,j) = sum(idx);
    Nhat(1,j) = sum(Nexit(1,1:j));
end

end
