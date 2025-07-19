nx = 25;
ny = 8;
nz = 25;

dataFolder = '/blue/bala1s/krishnap.kalivel/COVID/Sou2Snk_SFR/room_31/dist_1m/dt_060s/Data';
pFolder = '/blue/bala1s/krishnap.kalivel/COVID/Sou2Snk_SFR/room_31/dist_1m/dt_060s/temp';

%mkdir(pFolder);
%pFolder = '/blue/bala/krishnap.kalivel/COVID/room_31/Source2Sink';
ctr = 1;
for i=1:nx
    for j=1:ny
        for k=1:nz
            text = strcat('Source_nxi_',sprintf('%03d',round(i)),...
                '_nyi_',sprintf('%03d',round(j)),...
                '_nzi_',sprintf('%03d',round(k)),'.mat');
            filename = strcat(dataFolder,'/',text);
            
            if ~isfile(filename)
                text2 = strcat('folder_nxi_',sprintf('%03d',round(i)),...
                '_nyi_',sprintf('%03d',round(j)),...
                '_nzi_',sprintf('%03d',round(k)))
                folderloc = strcat(pFolder,'/',text2);
                if ctr
                    fid = fopen('mscript','wt');
                    ctr = 0;
                else
                    fid = fopen('mscript', 'a+');
                end
                fprintf(fid, 'cd ');
                fprintf(fid, folderloc);
                fprintf(fid, '\n');
                fprintf(fid, 'sbatch run.sh');
                fprintf(fid,'\n');
                fclose(fid);
            end
        end
    end
end

                
