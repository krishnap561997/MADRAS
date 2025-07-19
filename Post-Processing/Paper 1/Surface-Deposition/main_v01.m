clc, clear all, close all

% This routine extracts the location where particles deposit on room
% surface

%% Add path
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %s %s %s %f %f %f','Delimiter','\n');

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

%% Build the large matrix with information on unique_id, diameter, time of 
% deposit, surface of deposit, x,y,z coordinates of deposit
tic
disp('-----------------------------------')
disp('Building large matrix')
[id,dp,time,surf,loc] = deposit_v01(nfile);
toc
fprintf('\n')

%% Remove unaccounted particles that have exited the domain but do not appear in the parout files
tic
disp('-----------------------------------')
disp('Sorting')
[x_s,y_s,z_s,unique_diam,unique_surf,x,y,z] = sort_v01(dp,time,surf,loc);
toc
fprintf('\n')

%% Remove unaccounted particles that have exited the domain but do not appear in the parout files
tic
disp('-----------------------------------')
disp('Ploting')
plot_v01(x_s,y_s,z_s,unique_diam,unique_surf,x,y,z,fpath,fontsize);
toc
fprintf('\n')

tic
disp('-----------------------------------')
disp('saving file')
save(fullfile(fpathw, 'Deposit'))
disp('File saved')
toc

disp('-----------------------------------')
disp('Job terminated successfully')
a+b+c

