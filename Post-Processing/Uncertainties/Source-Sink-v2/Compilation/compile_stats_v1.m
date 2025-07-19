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

text = strcat(['CompiledData_Dist_',num2str(dist),'m.mat']);
load(fullfile(spath,text));
t_min = size(Nhat_cum{end},2);
t = 0:(t_min-1);

%% Total sigma
for i=1:numel(Nhat_cum)
    tmp = Nhat_cum{i};
    tmp2 = tmp(:,1:t_min);
    
    N_cum(:,i,:) = reshape(tmp2,[size(tmp2,1) 1 size(tmp2,2)]);
end

%% Process Data
N2_cum = reshape(mean(N_cum,2),[size(N_cum,1),size(N_cum,3)]) ;
N1_cum = reshape(mean(N_cum,1),[size(N_cum,2),size(N_cum,3)]) ;
N0_cum = reshape(N_cum, [size(N_cum,1).*size(N_cum,2),size(N_cum,3)]);

N1_std_v2 = reshape(std(N_cum,1,2),[size(N_cum,1),size(N_cum,3)]);
N2_std_v2 = reshape(std(N_cum,1,1),[size(N_cum,2),size(N_cum,3)]);
N_mean = mean(N0_cum,1);
N_std = std(N0_cum,1,1);
N2_mean = mean(N2_cum,1);
N2_std = std(N2_cum,1,1);
N1_std_calc = sqrt(N_std.^2 - N2_std.^2);
N1_std = std(N1_cum,1,1);
N1_std_v2 = mean(N1_std_v2,1);
N2_std_v2 = mean(N2_std_v2,1);

%% Save Data
fname = strcat(['StatsData_Dist_',num2str(dist),'m.mat']);
save(fullfile(spath,fname), 'N_mean', 'N_std','N2_std','N1_std','N1_std_calc','N1_std_v2','N2_std_v2');
