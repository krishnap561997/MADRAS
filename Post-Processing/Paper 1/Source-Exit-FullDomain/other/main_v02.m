clc, clear all, close all
%% Initial particles
% This routine computes the residence time of particles from a prescribed
% discrete ejection (such as a cough)

nfile = 10
%% Add path
% 'test1'
%addpath(genpath('C:\Users\nadim.zgheib\OneDrive - University of Florida\Covid_Framework\room_26'));
 addpath(genpath('/blue/bala/nzgheib/COVID/room_26/par'));
% 'test2'

%% Extract exit times for all tagged particles
% tic
[id,dp,time,x,y,z,exit_surf] = record_times_v02(nfile);
% toc
% 'test5'

%% Find the exit times for each particle size and interpolate for other
%  sizes in distribution
tic
disp('-----------------------------------')
disp('exit_time')
[Nwall,Nfloor,Noutlet,Ntot,unique_dp,unique_time] = exit_time_v02(dp,time,exit_surf);
toc

tic
disp('-----------------------------------')
disp('saving file')
save('NtotNexit_v03','Ntot','Nwall','Nfloor','Noutlet','unique_dp','unique_time')
disp('File saved')
toc
% plotting(Nexit,Ntot,unique_dp,unique_time)
