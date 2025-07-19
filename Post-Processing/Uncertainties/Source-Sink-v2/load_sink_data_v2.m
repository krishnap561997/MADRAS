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

function [id,dp,tm] = load_sink_data_v2(nxi,nyi,nzi,fpath)

name = strcat(['tm_nxi_',sprintf('%03d',round(nxi)),...
            '_nyi_',sprintf('%03d',round(nyi)),...
            '_nzi_',sprintf('%03d',round(nzi)),'.mat']);
load(fullfile(fpath,name));
tm = dum;

name = strcat(['dp_nxi_',sprintf('%03d',round(nxi)),...
            '_nyi_',sprintf('%03d',round(nyi)),...
            '_nzi_',sprintf('%03d',round(nzi)),'.mat']);
load(fullfile(fpath,name));
dp = dum;

name = strcat(['id_nxi_',sprintf('%03d',round(nxi)),...
            '_nyi_',sprintf('%03d',round(nyi)),...
            '_nzi_',sprintf('%03d',round(nzi)),'.mat']);
load(fullfile(fpath,name));
id = dum;

unique_dp = unique(dp,'sorted');
cond_dp = dp <= unique_dp(2);
tm = tm(cond_dp);
dp = dp(cond_dp);
id = id(cond_dp);

