clc, clear all, close all
dT1_max = 45;
dT2_max = 45;
d = 6;
T0_min = 0;
T0_max = 45;
ACH_min = 0.3;
ACH_max = 30;
parentFolder = '/blue/bala/krishnap.kalivel/COVID/HigherACHBetter/dist6m/temp';
for i=1:dT1_max
    for j=0:dT2_max
        dirtext = strcat(['folder_dT1_',sprintf('%03d',round(i)),...
                                '_dT2_',sprintf('%03d',round(j))])
        mkdir(parentFolder,dirtext);
        dir_path = strcat([parentFolder,'/',dirtext]);
        toto = char(strip(strip(dir_path,'left',''''),'right',''''));
        cd(toto);
        fid = fopen('input.txt','wt');
        fprintf(fid,"'/blue/bala/krishnap.kalivel/COVID/HigherACHBetter/dist6m/Data_dT1'");
        fprintf(fid,'\n');
        fprintf(fid,num2str(d));
        fprintf(fid,'\n');
        fprintf(fid,num2str(i));
        fprintf(fid,'\n');
        fprintf(fid,num2str(j));
        fprintf(fid,'\n');
        fprintf(fid,num2str(T0_min));
        fprintf(fid,'\n');
        fprintf(fid,num2str(T0_max));
        fprintf(fid,'\n');
        fprintf(fid,num2str(ACH_min));
        fprintf(fid,'\n');
        fprintf(fid,num2str(ACH_max));
        fprintf(fid,'\n');
        fclose(fid);

        copyfile('/home/krishnap.kalivel/MATLAB_src/Higher-ACH-Better/Main-Codes/v3/main_v3_dT1','main_v3')
        copyfile('/home/krishnap.kalivel/MATLAB_src/Higher-ACH-Better/Main-Codes/v3/run.sh','run.sh')
        cd('/blue/bala/krishnap.kalivel/COVID/HigherACHBetter/dist6m/temp');
        if(i==1 && j==0)
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

