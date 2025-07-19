clear all; close all; clc;
%location to store data, id_dp_tm location, file size, tlim, room
%dimensions, sphere diameter, distance between source and sink
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %s %f %f %f %f','Delimiter','\n');

fpath           = char(strip(strip(C{1},'left',''''),'right',''''))       %Save location
spath           = char(strip(strip(C{2},'left',''''),'right',''''))       %id_dp_tm location

nx              = C{3}  
ny              = C{4}     
nz              = C{5}      
dist            = C{6}
ctr=0;
for i=1:nx
    for j=1:ny
	tic
        text = strcat('FigureData_Dist_',num2str(dist),...
            'm_nx_', sprintf('%02d', round(i)),...
            '_ny_', sprintf('%02d', round(j)),'.mat');
        load(fullfile(spath, text));
        disp(text)
        if ctr==0
            Nexit_cum = Nexit_c;
            Nhat_cum = Nhat_c;
            sinks_loc = sinks_loc_c;
            source_loc = source_loc_c;
            ctr = size(sinks_loc,1);
        else
            for k = 1:numel(Nexit_cum)
                tmp = [Nexit_cum{k}; Nexit_c{k}];
                Nexit_cum{k} = tmp;
                tmp = [Nhat_cum{k}; Nhat_c{k}];
                Nhat_cum{k} = tmp;
                sinks_loc = [sinks_loc; sinks_loc_c];
                source_loc = [source_loc; source_loc_c];
            end
        end
        toc
    end 
end
text = strcat('CompiledData_Dist_',num2str(dist),'m.mat');
save(fullfile(spath,text),'Nhat_cum','Nexit_cum','sinks_loc','source_loc','-v7.3');
disp('Data Saved');

