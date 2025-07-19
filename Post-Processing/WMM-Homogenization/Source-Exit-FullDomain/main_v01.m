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

clc, clear all, close all
%% Initial particles
% This routine computes the residence time of particles for entire room

%% Add path
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %s %s %s %f %f %f %f %f %f','Delimiter','\n');

par     = char(strip(strip(C{1},'left',''''),'right',''''))
parout  = char(strip(strip(C{2},'left',''''),'right',''''))
fpath   = char(strip(strip(C{3},'left',''''),'right',''''))
fpathw  = char(strip(strip(C{4},'left',''''),'right',''''))
%parent  = char(strip(strip(C{5},'left',''''),'right',''''));
fclose(fileID);

addpath(genpath(par))
addpath(genpath(parout))
%addpath(genpath(parent))
nfile   = C{5}
ACH     = C{6}
fontsize= C{7}
rlenx   = C{8}
rleny   = C{9}
rlenz   = C{10}

%% Extract exit times for all tagged particles
tic
disp('-----------------------------------')
disp('Recording exit times and surfaces')
[id,dp,time,surf,loc] = record_times_v01(nfile);
toc
fprintf('\n')

%% Remove unaccounted particles that have exited the domain but do not appear in the parout files
tic
disp('-----------------------------------')
disp('verify_remove')
[id,dp,time,surf,loc] = verify_remove_01(id,dp,time,surf,loc,nfile);
toc
fprintf('\n')

tic
disp('-----------------------------------')
disp('exit_time')
[Nwall,Nfloor,Noutlet,Nceil,Nexit,Ntot,unique_dp] = exit_time_v01(dp,time,nfile,surf,loc,rlenx,rlenz);
toc
fprintf('\n')

%% Compute settling velocity
rhop = 1e3; rg = 9.8;   rmu = 1.81e-5;  
rhoa = 1.2041;   rhoh = rhop/rhoa;
beta = 3/(2*rhoh+1);    nu = rmu/rhoa;
Vs   = 0*unique_dp;  
for i=1:10
    rep     = Vs .* unique_dp ./ nu;
    tau     = unique_dp.^2*(rhoh+0.5)./(18*rmu*(1+0.15*rep.^0.687));
    Vs      = tau.*(1-beta)*rg;
end
H = 3.2;

%% Plot
tic
disp('-----------------------------------')
disp('plotting')
lambda = plot_v01(Nfloor,Noutlet,Nwall,Nceil,Nexit,Ntot,unique_dp,fpath,Vs,ACH,H,fontsize);
toc
fprintf('\n')

tic
disp('-----------------------------------')
disp('saving file')
%save('NtotNexit_v01')
save(fullfile(fpathw, 'NtotNexit_v01'))
disp('File saved')
toc

disp('-----------------------------------')
disp('Job terminated successfully')
a+b+c
