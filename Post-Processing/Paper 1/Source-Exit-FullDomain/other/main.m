clc, clear all, close all
%% Initial particles
% This routine computes the residence time of particles from a prescribed
% discrete ejection (such as a cough)

nfile = 281
%% Add path
% 'test1'
% addpath(genpath('C:\Users\nadim.zgheib\OneDrive - University of Florida\Covid_Framework\room_25'));
%addpath(genpath('C:\Users\nadim.zgheib\OneDrive - University of Florida\Covid_Framework\room_26'));
 addpath(genpath('/blue/bala/nzgheib/COVID/room_26/par'));
% 'test2'

%% Extract exit times for all tagged particles
% tic
[id,dp,time] = record_times(nfile);
% toc
% 'test5'

%% Find the exit times for each particle size and interpolate for other
%  sizes in distribution
[Nexit,Ntot,unique_dp,unique_time] = exit_time(dp,time);

save('NtotNexit','Ntot','Nexit','unique_dp','unique_time')
disp('File saved')
