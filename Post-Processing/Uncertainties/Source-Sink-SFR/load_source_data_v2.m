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

function [id,dp,tm,x0,y0,z0] = load_source_data_v2(nxi,nyi,nzi,fpath)

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

% idx = (tm>nfile-time_limit);
% tm(idx) = [];
% dp(idx) = [];
% id(idx) = [];

% Considering only 0.2 and 1 micron particles
unique_dp = unique(dp,'sorted');
cond_dp = dp <= unique_dp(2);
tm = tm(cond_dp);
dp = dp(cond_dp);
id = id(cond_dp);

% Loading x,y,z locations if needed
name = strcat(['loc_nxi_',sprintf('%03d',round(nxi)),...
            '_nyi_',sprintf('%03d',round(nyi)),...
            '_nzi_',sprintf('%03d',round(nzi)),'.mat']);
if isfile(fullfile(fpath,name))
    load(fullfile(fpath,name));
    x0 = dum(:,1);
    y0 = dum(:,2);
    z0 = dum(:,3);
    
    x0(idx) = [];
    y0(idx) = [];
    z0(idx) = [];
    
    x0 = x0(cond_dp);
    y0 = y0(cond_dp);
    z0 = z0(cond_dp);
else
    x0 = [];
    y0 = [];
    z0 = [];
end



