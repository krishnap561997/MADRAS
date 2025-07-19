clear; clc; close all;
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %s %f %f %f %f','Delimiter','\n');

fpath           = char(strip(strip(C{1},'left',''''),'right',''''))       %Save location
spath           = char(strip(strip(C{2},'left',''''),'right',''''))       %id_dp_tm location

nx              = C{3}
ny              = C{4}
nz              = C{5}
dist            = C{6}


ctr = 1;
for i=1:nx
    for j=1:ny
        tic
        for k=1:nz
            text = strcat(['Source_nxi_',sprintf('%03d',round(i)),...
                '_nyi_',sprintf('%03d',round(j)),...
                '_nzi_',sprintf('%03d',round(k))])
	    if isfile(fullfile(fpath,text))
		    disp('No source to sink combinations for this source');
		    continue;
	    end
            load(fullfile(fpath, text));
            source_loc = [i,j,k].*ones(size(sinks_loc));
            if ctr == 1
                source_loc_c = source_loc;
                sinks_loc_c = sinks_loc;
            else
                source_loc_c = [source_loc_c; source_loc];
                sinks_loc_c = [sinks_loc_c; sinks_loc];
            end

            
            for l1=1:size(Nhat_cum,1)
                for l2 = 1:size(Nhat_cum,2)
                    tmp = Nhat_cum{l1,l2};
                    Nhat_ss_cum(ctr,l2) = tmp(end)./Ntot_cum(l1,l2);
                    time_max(ctr,l2) = numel(tmp)-1;
                end
		ctr = ctr+1;
            end

        end
        toc
    end
end
text = strcat('SteadyStateData_Dist_',num2str(dist),'m.mat');
save(fullfile(spath,text),'Nhat_ss_cum','sinks_loc_c','source_loc_c','time_max');
disp('Data Saved');
        
