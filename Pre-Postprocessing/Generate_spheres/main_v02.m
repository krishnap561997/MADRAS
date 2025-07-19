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
% The routine computes residence time of particles from a finite source

%% Add path
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %s %f %f %f %f %f %f %f','Delimiter','\n');

par     = char(strip(strip(C{1},'left',''''),'right',''''))
fpath   = char(strip(strip(C{2},'left',''''),'right',''''))
%parent  = char(strip(strip(C{5},'left',''''),'right',''''));
fclose(fileID);

addpath(genpath(par))
%addpath(genpath(parent))
nfile           = C{3}
sphere_diameter = C{4}
nxi             = C{5}
nzi             = C{6}
rlenx           = C{7}
rleny           = C{8}
rlenz           = C{9}

%% Read parameters
[xsphere,ysphere,zsphere] = center_location(nxi,nzi,sphere_diameter,rlenx,rleny,rlenz)
%% Identify particles within the sphere of interest for all output files
tic
disp('-----------------------------------')
disp('Building large vectors')
[id,dp,tm] = sph_tag(sphere_diameter,xsphere,ysphere,zsphere,nfile);
toc
fprintf('\n')

%% Save vectors in separate files
tic
disp('-----------------------------------')
disp('Saving the workspace')
for nyi=1:numel(ysphere)
    tic
    disp('-----------------------------------')
    text = strcat(['nxi=',sprintf('%03d',round(nxi)),'nyi=',sprintf('%03d',round(nyi))...
        ,'_nzi=',sprintf('%03d',round(nzi)),'_sphere_data']);
    disp(text)
    %% Saving the particle id
    text = strcat(['id_nxi_',sprintf('%03d',round(nxi)),'_nyi_',sprintf('%03d',round(nyi)),...
        '_nzi_',sprintf('%03d',round(nzi))]);
    dum = id{nyi};
    save(fullfile(fpath, text),'dum');
    disp('id Saved')
    
    %% Saving the particle diam
    text = strcat(['dp_nxi_',sprintf('%03d',round(nxi)),'_nyi_',sprintf('%03d',round(nyi)),...
        '_nzi_',sprintf('%03d',round(nzi))]);
    dum = dp{nyi};
    save(fullfile(fpath, text),'dum');
    disp('diameter Saved')
    
    %% Saving the particle time
    text = strcat(['tm_nxi_',sprintf('%03d',round(nxi)),'_nyi_',sprintf('%03d',round(nyi)),...
        '_nzi_',sprintf('%03d',round(nzi))]);
    dum = tm{nyi};
    save(fullfile(fpath, text),'dum');
    disp('time Saved')
    toc

    tic
    disp('-----------------------------------')
    disp('------Making the data unique-------')
    [id_unique,dp_unique,tm_unique] = make_unique(id{nyi},dp{nyi},tm{nyi});
    text = strcat(['Unique: nxi=',sprintf('%03d',round(nxi)),'nyi=',sprintf('%03d',round(nyi))...
        ,'_nzi=',sprintf('%03d',round(nzi)),'_sphere_data']);
    disp(text)
    %% Saving the particle id
    text = strcat(['unq_id_nxi_',sprintf('%03d',round(nxi)),'_nyi_',sprintf('%03d',round(nyi)),...
        '_nzi_',sprintf('%03d',round(nzi))]);
    dum = id_unique;
    save(fullfile(fpath, text),'dum');
    disp('id Saved')
    
    %% Saving the particle diam
    text = strcat(['unq_dp_nxi_',sprintf('%03d',round(nxi)),'_nyi_',sprintf('%03d',round(nyi)),...
        '_nzi_',sprintf('%03d',round(nzi))]);
    dum = dp_unique;
    save(fullfile(fpath, text),'dum');
    disp('diameter Saved')
    
    %% Saving the particle time
    text = strcat(['unq_tm_nxi_',sprintf('%03d',round(nxi)),'_nyi_',sprintf('%03d',round(nyi)),...
        '_nzi_',sprintf('%03d',round(nzi))]);
    dum = tm_unique;
    save(fullfile(fpath, text),'dum');
    disp('time Saved')
    toc

end
disp('id, dp, & tm Saved')
toc

disp('-----------------------------------')
disp('Job terminated successfully')
a+b+c

