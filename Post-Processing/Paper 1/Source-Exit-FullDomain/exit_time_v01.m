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

function [Nwall,Nfloor,Noutlet,Nceil,Nexit,Ntot,unique_dp] = exit_time_v01(dp,time,nfile,surf,loc,rlenx,rlenz)

%% Compute the number of particle sizes involved
unique_dp   = unique(dp);


Ntot    = zeros(numel(unique_dp),1);
Nwall   = zeros(numel(unique_dp),nfile);
Nfloor  = zeros(numel(unique_dp),nfile);
Noutlet = zeros(numel(unique_dp),nfile);
Nceil   = zeros(numel(unique_dp),nfile);
Nexit   = zeros(numel(unique_dp),nfile);

for i=1:numel(unique_dp) % Number of particle sizes
    idx         = (dp==unique_dp(i));
    Ntot(i,1)   = sum(idx);   
    for j = 1:nfile % Number of time instances
        idx         = (dp==unique_dp(i) & (time<=j-1) & surf==3 & (loc(:,1)-rlenx/2).^2 + (loc(:,3)-rlenz/2).^2>0.8^2);
        Nceil(i,j)  = sum(idx);       
        
        idx         = (dp==unique_dp(i) & (time<=j-1) & surf==1);
        Nfloor(i,j)  = sum(idx);       
        
        idx         = (dp==unique_dp(i) & (time<=j-1) & (surf==7 | (surf==3 & (loc(:,1)-rlenx/2).^2 + (loc(:,3)-rlenz/2).^2<=0.8^2)));
        Noutlet(i,j)  = sum(idx);       
        
        idx         = (dp==unique_dp(i) & (time<=j-1) & surf>1 & surf<7 & surf~=3);
        Nwall(i,j)  = sum(idx);       

        Nexit(i,j)  = Nceil(i,j) + Nfloor(i,j) + Noutlet(i,j) + Nwall(i,j);
    end
end

