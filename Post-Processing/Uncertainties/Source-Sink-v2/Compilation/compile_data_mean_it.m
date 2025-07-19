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

for i=1:nx
    
    for j=1:ny
        ctr = 1;
        for k=1:nz
            tic
            text = strcat(['Source_nxi_',sprintf('%03d',round(i)),...
                '_nyi_',sprintf('%03d',round(j)),...
                '_nzi_',sprintf('%03d',round(k))]);
            load(fullfile(fpath, text));
	    if numel(sinks_loc)==0
		disp('Skipped due to no Src-Snk combinations');
	        continue;
	    end
            disp(text)
            source_loc = [i,j,k].*ones(size(sinks_loc));
            for l1=1:size(sinks_loc,1)
                sinks_loc_c(ctr,:) = sinks_loc(l1,:);
                source_loc_c(ctr,:) = [i,j,k];
                if ctr ==1
                    Nexit_c = cell(size(Nexit_cum,2),1);
                    Nhat_c = cell(size(Nexit_cum,2),1);
                    for l2 = 1:numel(Nexit_c)
                        Nexit_c{l2} = Nexit_cum{l1,l2}./Ntot_cum(l1,l2);
                        Nhat_c{l2} = Nhat_cum{l1,l2}./Ntot_cum(l1,l2);
                    end
                else

                    for l2 = 1:numel(Nexit_c)
                        tmp1 = Nexit_c{l2};
                        tmp1(ctr,:) = Nexit_cum{l1,l2}./Ntot_cum(l1,l2);

                        tmp2 = Nhat_c{l2};
                        tmp2(ctr,:) = Nhat_cum{l1,l2}./Ntot_cum(l1,l2);
                        Nexit_c{l2} = tmp1;
                        Nhat_c{l2} = tmp2;
                    end
                end

                ctr = ctr+1; 
            end
            toc
        end
        text = strcat('FigureData_Dist_',num2str(dist),...
            'm_nx_', sprintf('%02d', round(i)),...
            '_ny_', sprintf('%02d', round(j)),'.mat');
        save(fullfile(spath,text),'Nhat_c','Nexit_c','sinks_loc_c','source_loc_c');
	disp('Data Saved');
        clear tmp1 tmp2 Nhat_c Nexit_c sinks_loc_c source_loc_c;
    end 
end
