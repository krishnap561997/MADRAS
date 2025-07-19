clear all; close all; clc;
%location to store data, id_dp_tm location, file size, tlim, room
%dimensions, sphere diameter, distance between source and sink
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %s %f %f %f %f','Delimiter','\n');

fpath           = char(strip(strip(C{2},'left',''''),'right',''''));        %Save location
spath           = pwd;

nx              = C{3}  
ny              = C{4}     
nz              = C{5}      
dist            = C{6}
ctr  = 0;

for k=1:ceil(nz/2)  
    for i=k:ceil(nx/2)
        for j=1:ny  
	    
            text = strcat(['Source_8f_nxi_',sprintf('%03d',round(i)),...
                '_nyi_',sprintf('%03d',round(j)),...
                '_nzi_',sprintf('%03d',round(k)),'.mat']);

	    try
                load(fullfile(fpath, text));
                disp(text)
	    catch ME
		disp(text)
		disp('Skipped');
		continue;
   	    end
            
            if i==1 && j==1 && k==1
                Nsfr_cum = Nsfr;
                src_loc_cum = src_loc_c;
                snk_loc_cum = snk_loc_c;
            else
                Nsfr_cum = [Nsfr_cum; Nsfr];
                src_loc_cum = [src_loc_cum; src_loc_c];
                snk_loc_cum = [snk_loc_cum; snk_loc_c];
            end	    
        end
    end
end

tmp = std(Nsfr_cum,0,3);
tmp(isnan(tmp)) = 0;
sigma1 = mean(tmp,1);

tmp = skewness(Nsfr_cum,0,3);
tmp(isnan(tmp)) = 0;
skew1 = mean(tmp,1);

tmp = kurtosis(Nsfr_cum,0,3);
tmp(isnan(tmp)) = 0;
kurt1 = mean(tmp ,1);


tmp = std(Nsfr_cum,0,1);
tmp(isnan(tmp)) = 0;
sigma2 = mean(tmp,3);

tmp = skewness(Nsfr_cum,0,1);
tmp(isnan(tmp)) = 0;
skew2 = mean(tmp,3);

tmp = kurtosis(Nsfr_cum,0,1);
tmp(isnan(tmp)) = 0;
kurt2 = mean(tmp ,3);

avg = mean(mean(Nsfr_cum,1),3);


text = strcat(['Stats_Dist_',num2str(dist),'m.mat']);
save(fullfile(spath,text),'avg','sigma1','sigma2','kurt1','kurt2','skew1','skew2');

%text = strcat(['CompiledData_Dist_',num2str(dist),'m.mat']);
%save(fullfile(spath,text),'Nsfr_cum','snk_loc_cum','src_loc_cum');

