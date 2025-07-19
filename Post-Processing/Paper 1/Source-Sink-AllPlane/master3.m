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
time_limit = 1600;
nxt = 25;    nyt = 8;     nzt = 25;
rlenx = 10; rleny = 3.2; rlenz = 10;
sphere_diameter = 0.4; 
nxi = 13;
nyi = 3;
nzi = 13;
parentFolder = '/blue/bala/krishnap.kalivel/COVID/Sou2Snk_AllPlane/room_30';
text1 = strcat(['Source_nxi_',sprintf('%03d',round(nxi)),...
                                '_nyi_',sprintf('%03d',round(nyi)),...;
                                '_nzi_',sprintf('%03d',round(nzi))]);
mkdir(parentFolder,text1);
parentFolder = strcat(parentFolder,'/',text1);
X = find_symm(nxi, nzi, rlenx, rlenz, sphere_diameter);

for i=1:size(X,1)
    nxi = X(i,1);
    nzi = X(i,2);
    dirtext = strcat(['Source_nxi_',sprintf('%03d',round(nxi)),...
                                '_nyi_',sprintf('%03d',round(nyi)),...;
                                '_nzi_',sprintf('%03d',round(nzi))]);
    mkdir(parentFolder,dirtext)
    dir_path = strcat([parentFolder,'/',dirtext]);
    toto = char(strip(strip(dir_path,'left',''''),'right',''''));
    cd(toto)
    parentFolder2 = strcat(parentFolder, '/', dirtext);
    for i1 = 1:nxt
        for j=1:nyt
            for k=1:nzt
            dirtext2 = strcat(['Sink_', sprintf('%03d',round(i1)),...
                '_',sprintf('%03d',round(j)),'_',sprintf('%03d',round(k))]);
            mkdir(parentFolder2,dirtext2)
            dir_path = strcat([parentFolder2,'/',dirtext2]);
            toto = char(strip(strip(dir_path,'left',''''),'right',''''));
            cd(toto)
            text2 = strcat("'",parentFolder2, "'");
            fid = fopen('input.txt','wt');
            fprintf(fid,text2);
            fprintf(fid,'\n');
            fprintf(fid,"'/blue/bala/krishnap.kalivel/COVID/room_30/id_dp_tm'");
            fprintf(fid,'\n');
            fprintf(fid,num2str(nfile));
            fprintf(fid,'\n');
            fprintf(fid,num2str(time_limit));
            fprintf(fid,'\n');
            fprintf(fid,num2str(nxi));
            fprintf(fid,'\n');
            fprintf(fid,num2str(nyi));
            fprintf(fid,'\n');
            fprintf(fid,num2str(nzi));
            fprintf(fid,'\n');
            fprintf(fid,num2str(rlenx));
            fprintf(fid,'\n');
            fprintf(fid,num2str(rleny));
            fprintf(fid,'\n');
            fprintf(fid,num2str(rlenz));
            fprintf(fid,'\n');
            fprintf(fid,num2str(sphere_diameter));
            fprintf(fid,'\n');
            fprintf(fid,num2str(i1));
            fprintf(fid,'\n');
            fprintf(fid,num2str(j));
            fprintf(fid,'\n');
            fprintf(fid,num2str(k));
            fprintf(fid,'\n');
            fclose(fid);

            tmp = '/home/krishnap.kalivel/MATLAB_src/Post-Processing/Source-Sink-AllPlane';
            text = strcat(tmp,'/main_sou2snk');
            copyfile(text,'main_sou2snk')
            text = strcat(tmp,'/run.sh');
            copyfile(text,'run.sh')
            cd(parentFolder2);
            if(i1==1 && j ==1 && k==1)
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

end


