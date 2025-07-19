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

function [lambda] = plot_v01(Nfloor,Noutlet,Nwall,Nceil,Nexit,Ntot,unique_dp,fpath,Vs,ACH,H,fontsize)

for i=1:numel(unique_dp)
    CM = jet(i);
end
time_limit = numel(Nfloor(1,:))-1;

figure_v01(CM,unique_dp,0:time_limit,Nexit./Ntot,'t\,(s)'...
    ,'C_{exit}','Nexit',fpath,fontsize)
figure_v01(CM,unique_dp,0:time_limit,Nfloor./Ntot,'t\,(s)'...
    ,'C_{floor}','Nfloor',fpath,fontsize)
figure_v01(CM,unique_dp,0:time_limit,Noutlet./Ntot,'t\,(s)'...
    ,'C_{outlet}','Noutlet',fpath,fontsize)
figure_v01(CM,unique_dp,0:time_limit,Nwall./Ntot,'t\,(s)'...
    ,'C_{walls}','Nwall',fpath,fontsize)
figure_v01(CM,unique_dp,0:time_limit,Nceil./Ntot,'t\,(s)'...
    ,'C_{ceiling}','Nceil',fpath,fontsize)

figure_v01(CM,unique_dp,0:time_limit,Nfloor./Nexit,'t\,(s)'...
    ,'C_{floor}/C_{exit}','Nfloor_by_Nexit',fpath,fontsize)
figure_v01(CM,unique_dp,0:time_limit,Noutlet./Nexit,'t\,(s)'...
    ,'C_{outlet}/C_{exit}','Noutlet_by_Nexit',fpath,fontsize)
figure_v01(CM,unique_dp,0:time_limit,Nwall./Nexit,'t\,(s)'...
    ,'C_{walls}/C_{exit}','Nwall_by_Nexit',fpath,fontsize)
figure_v01(CM,unique_dp,0:time_limit,Nceil./Nexit,'t\,(s)'...
    ,'C_{ceiling}/C_{exit}','Nceil_by_Nexit',fpath,fontsize)

%% Exponential fitting
[lambda] = exp_fitting_v01(0:time_limit,Nexit./Ntot,...
    't','1-C/C_0',fpath,unique_dp,Vs,ACH,H,fontsize);
