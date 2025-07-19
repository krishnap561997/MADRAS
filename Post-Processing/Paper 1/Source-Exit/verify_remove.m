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

function [id,dp,time,time2,exit_surf] = verify_remove(id,dp,time,time2,exit_surf,nfile)
infile = strcat(sprintf('par%05d.h5\n', nfile));
%disp('-----------------------------------')
disp(['Reading ',infile])
time_nf        = h5read(infile,'/g1/time');% single
data_id_nf     = h5read(infile,'/g1/tag'); % unit32
data_dp_nf     = h5read(infile,'/g1/dp');  % single
data_loc_nf    = h5read(infile,'/g1/loc'); % single

[Lia2,~] = ismember(id,data_id_nf);

toto = (time2).*(Lia2);
if (sum(toto)~=0)
    disp('Warning: some particles that exited domain are still airborne')
end

toto2 = time2+Lia2;
idx2 = find(toto2==0);
if (sum(numel(idx2))~=0)
    disp([num2str(sum(numel(idx2))),' unaccounted particles out of ',num2str(sum(numel(time2))),' have been removed'])
end

id(idx2)        = [];
dp(idx2)        = [];
time(idx2)      = [];
time2(idx2)     = [];
exit_surf(idx2) = [];

