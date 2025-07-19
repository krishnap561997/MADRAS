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

% Only source arrays are unique

function [src2snkDt] = src2snk_v02(src_id,src_tm,snk_id,snk_tm)

src2snkDt = 0*snk_tm + 1e9;
[tf,idx] = ismember(snk_id,src_id);
src2snkDt(tf) = snk_tm(tf) - src_tm(idx(idx ~=0));

for i=1:size(src2snkDt,1)
    if src2snkDt(i)<0
        src2snkDt(i)=1e9;
    end
end


