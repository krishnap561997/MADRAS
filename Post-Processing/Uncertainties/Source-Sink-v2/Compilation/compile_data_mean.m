clear all; close all; clc;
%location to store data, id_dp_tm location, file size, tlim, room
%dimensions, sphere diameter, distance between source and sink
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %s %f %f %f %f','Delimiter','\n');

fpath           = char(strip(strip(C{1},'left',''''),'right',''''));        %Save location
spath           = char(strip(strip(C{2},'left',''''),'right',''''));      %id_dp_tm location

nx              = C{3}  
ny              = C{4}     
nz              = C{5}      
dist            = C{6}
ctr  = 0;
for i=1:nx
    for j=1:ny
        for k=1:nz
	    tic
            text = strcat(['Source_nxi_',sprintf('%03d',round(i)),...
                '_nyi_',sprintf('%03d',round(j)),...
                '_nzi_',sprintf('%03d',round(k))]);
            load(fullfile(fpath, text));
            disp(text)
            for l1=1:size(sinks_loc,1)
                if ctr==0                
                        Nexit = cell(size(Nexit_cum,2),1);
                        Nhat = cell(size(Nexit_cum,2),1);
			Nhat_all = cell(size(Nexit_cum,2),1);
                end
                for l2 = 1:numel(Nexit)
                    if ctr==0 
                        Nexit{l2} = Nexit_cum{l1,l2}./Ntot_cum(l1,l2);
                        Nhat{l2} = Nhat_cum{l1,l2}./Ntot_cum(l1,l2);
			Nhat_all{l2} = Nhat{l2};
                    else
                        Nexit{l2} = Nexit{l2} + Nexit_cum{l1,l2}./Ntot_cum(l1,l2);
                        Nhat{l2} = Nhat{l2} + Nhat_cum{l1,l2}./Ntot_cum(l1,l2);
			tmp = Nhat_all{l2};
                        tmp(ctr+1,:) = Nhat_cum{l1,l2}./Ntot_cum(l1,l2);
                        Nhat_all{l2} = tmp;
                    end
                end
                
                ctr = ctr+1;
            end
	    toc
        end
    end
end

for i=1:numel(Nexit)
    Nexit_avg{i} = Nexit{i}./ctr;
    Nhat_avg{i} = Nhat{i}./ctr;
end

text = strcat(['CompiledData_Dist_',num2str(dist),'m.mat']);
save(fullfile(spath,text),'Nhat_avg','Nhat_all','Nexit_avg','ctr');
