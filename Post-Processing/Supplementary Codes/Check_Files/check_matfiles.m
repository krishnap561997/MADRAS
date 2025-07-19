nsrc2snk = 10000;

dataFolder = '/blue/bala1s/krishnap.kalivel/COVID/Sou2Snk_WMM/room_34/AllSizes/Data';
pFolder = '/blue/bala1s/krishnap.kalivel/COVID/Sou2Snk_WMM/room_34/AllSizes/temp';

%mkdir(pFolder);
%pFolder = '/blue/bala/krishnap.kalivel/COVID/room_31/Source2Sink';
ctr = 1;
for i=1:nsrc2snk

            text = strcat(['Source_s2sno_',sprintf('%05d',round(i))]);
            % filename = strcat(dataFolder,'/',text);
            disp(num2str(i))
            try
                load(fullfile(dataFolder,text));
            catch ME
                text2 = strcat(['folder_',sprintf('%05d',round(i))])
                folderloc = fullfile(pFolder,text2);
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

                
