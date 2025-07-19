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
nfile = 2145;
time_limit_min = 1144;
time_gap = 200;
parentFolder ='/blue/bala/krishnap.kalivel/RANS/RANS_LGRoom/Sou2Snk/dist_2m/temp';
nxt = 19;    nyt = 5;     nzt = 29;
rlenx = 7.6; rleny = 2.3; rlenz = 11.8;
sphere_diameter = 0.4; dist = 2; tolerance = 0;
for i=1:nxt
    for j=1:nyt
        for k=1:nzt
            dirtext = strcat(['folder_nxi_',sprintf('%03d',round(i)),...
                                    '_nyi_',sprintf('%03d',round(j)),...;
                                    '_nzi_',sprintf('%03d',round(k))]);
            mkdir(parentFolder,dirtext)
            dir_path = strcat([parentFolder,'/',dirtext]);
            toto = char(strip(strip(dir_path,'left',''''),'right',''''));
            cd(toto)
            fid = fopen('input.txt','wt');
            fprintf(fid,"'/blue/bala/krishnap.kalivel/RANS/RANS_LGRoom/Sou2Snk/dist_2m/Data'");
            fprintf(fid,'\n');
            fprintf(fid,"'/blue/bala/krishnap.kalivel/RANS/RANS_LGRoom/id_dp_tm'");
            fprintf(fid,'\n');
            fprintf(fid,num2str(nfile));
            fprintf(fid,'\n');
            fprintf(fid,num2str(time_limit_min));
            fprintf(fid,'\n');
            fprintf(fid,num2str(time_gap));
            fprintf(fid,'\n');
            fprintf(fid,num2str(i));
            fprintf(fid,'\n');
            fprintf(fid,num2str(j));
            fprintf(fid,'\n');
            fprintf(fid,num2str(k));
            fprintf(fid,'\n');
            fprintf(fid,num2str(rlenx));
            fprintf(fid,'\n');
            fprintf(fid,num2str(rleny));
            fprintf(fid,'\n');
            fprintf(fid,num2str(rlenz));
            fprintf(fid,'\n');
            fprintf(fid,num2str(sphere_diameter));
            fprintf(fid,'\n');
            fprintf(fid,num2str(dist));
            fprintf(fid,'\n');
            fprintf(fid,num2str(tolerance));
            fclose(fid);

            copyfile('/home/krishnap.kalivel/MATLAB_src/Aerosol-Paper/Source-Sink-v2/main_sou2snk_v2','main_sou2snk_v2')
            copyfile('/home/krishnap.kalivel/MATLAB_src/Aerosol-Paper/Source-Sink-v2/run.sh','run.sh')
            cd('/blue/bala/krishnap.kalivel/RANS/RANS_LGRoom/Sou2Snk/dist_2m/temp');
            if(i==1 && j==1 &&  k==1)
               fid = fopen('mscript','wt');
            else
               fid = fopen('mscript', 'a+');
            end
            fprintf(fid,'cd ');
            fprintf(fid,toto);
            fprintf(fid,'\n');
            fprintf(fid,'sbatch run.sh');
            fprintf(fid,'\n');
            fclose(fid);
        end
    end
end



