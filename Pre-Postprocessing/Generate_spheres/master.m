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
nfile = 1813;
sphere_diameter = 0.2;
rlenx = 10;
rleny = 3.2;
rlenz = 10;
parentFolder = '/blue/bala/krishnap.kalivel/COVID/room_30/temp';
nxt = 25;    nzt = 25;
for i=1:nxt
    for k=1:nzt
        dirtext = strcat(['folder_nxi_',sprintf('%03d',round(i)),...
            '_nzi_',sprintf('%03d',round(k))]);
        mkdir(parentFolder,dirtext)
        dir_path = strcat([parentFolder,'/',dirtext]);
        toto = char(strip(strip(dir_path,'left',''''),'right',''''));
        cd(toto)
        fid = fopen('input.txt','wt');
        fprintf(fid,"'/blue/bala/nzgheib/COVID/room_30/par'");
        fprintf(fid,'\n');
        fprintf(fid,"'/blue/bala/krishnap.kalivel/COVID/room_30/id_dp_tm_2'");
        fprintf(fid,'\n');
        fprintf(fid,num2str(nfile));
        fprintf(fid,'\n');
        fprintf(fid,num2str(sphere_diameter));
        fprintf(fid,'\n');
        fprintf(fid,num2str(i));
        fprintf(fid,'\n');
        fprintf(fid,num2str(k));
        fprintf(fid,'\n');
        fprintf(fid,num2str(rlenx));
        fprintf(fid,'\n');
        fprintf(fid,num2str(rleny));
        fprintf(fid,'\n');
        fprintf(fid,num2str(rlenz));
        fclose(fid);
   
        copyfile('/home/krishnap.kalivel/MATLAB_src/Post-Processing/Generate_spheres/main_v01','main_v01'); 
        copyfile('/home/krishnap.kalivel/MATLAB_src/Post-Processing/Generate_spheres/run.sh','run.sh');
        cd('/blue/bala/krishnap.kalivel/COVID/room_30/temp');
        if(i==1 && k==1)
           fid = fopen('mscript','wt');
        else
           fid = fopen('mscript', 'a+');
        end
        fprintf(fid,'cd ')
        fprintf(fid,toto);
        fprintf(fid,'\n');
        fprintf(fid,'sbatch run.sh');
        fprintf(fid,'\n');
        fclose(fid);
    end
end

