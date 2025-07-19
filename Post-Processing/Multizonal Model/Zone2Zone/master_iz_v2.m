clear; clc; close all;

nfile = 4479;

nsrc2snk = 5000;

time_gap = 50;
parentFolder ='/blue/bala1s/krishnap.kalivel/COVID/IZModel/Code/Zone2Zone/temp';
rlenx = 10; rleny = 3.2; rlenz = 10;
outlet_x = 5; outlet_z = 5;
dw = 0.1;
t_start = 0:time_gap:1000;
ctr = 1;
for i=1:4
    for j=1:numel(t_start)
        dirtext = strcat(['folder_zone_',num2str(i),'_t_',sprintf('%04d',round(t_start(j)))]);
        mkdir(parentFolder,dirtext)
    
        dir_path = strcat([parentFolder,'/',dirtext]);
        toto = char(strip(strip(dir_path,'left',''''),'right',''''));
        cd(toto)
        fid = fopen('input.txt','wt');
        fprintf(fid,"'/orange/bala1s/krishnap.kalivel/COVID/room_34/h5/par'");
        fprintf(fid,'\n');
        %fprintf(fid,"'/blue/bala1s/rupal.patel/Results/room_one_HT6/data/id_dp_tm'");
        fprintf(fid,"'/blue/bala1s/krishnap.kalivel/COVID/IZModel/room_34/Data'");
        fprintf(fid,'\n');
        fprintf(fid,num2str(nfile));
        fprintf(fid,'\n');
        fprintf(fid,num2str(rlenx));
        fprintf(fid,'\n');
        fprintf(fid,num2str(rleny));
        fprintf(fid,'\n');
        fprintf(fid,num2str(rlenz));
        fprintf(fid,'\n');
        fprintf(fid,num2str(outlet_x));
        fprintf(fid,'\n');
        fprintf(fid,num2str(outlet_z));
        fprintf(fid,'\n');
        fprintf(fid,num2str(t_start(j)));
        fprintf(fid,'\n');
        fprintf(fid,num2str(i));
        fprintf(fid,'\n');
        fprintf(fid,num2str(dw));
        fprintf(fid,'\n');

        fclose(fid);
    
        copyfile('/blue/bala1s/krishnap.kalivel/COVID/IZModel/Code/Zone2Zone/main_zone2zone_v2','main_zone2zone_v2')
        copyfile('/blue/bala1s/krishnap.kalivel/COVID/IZModel/Code/Zone2Zone/run.sh','run.sh')
        cd('/blue/bala1s/krishnap.kalivel/COVID/IZModel/Code/Zone2Zone/temp');
        if(ctr==1)
           fid = fopen('mscript','wt');
           ctr = 0;
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
