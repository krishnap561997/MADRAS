%% Add Path and Load Input File. Contains par, parout file locations (NOT NEEDED FOR THIS),
clear; close; clc;
%location to store data, id_dp_tm location, file size, tlim, room
%dimensions, sphere diameter, distance between source and sink
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %s %f %f %f %f %f %f %f %f %f %f %f','Delimiter','\n');

fpath   = char(strip(strip(C{1},'left',''''),'right',''''));        %Save location
dpath     = char(strip(strip(C{2},'left',''''),'right',''''));      %id_dp_tm location

nfile           = C{3}      %Number of par files
time_gap        = C{4}      %Gap between ejection periods
nxi             = C{5}      %Source x,y,z indices
nyi             = C{6}
nzi             = C{7}
rlenx           = C{8}      %Room size, x, y, z
rleny           = C{9}
rlenz           = C{10}
sphere_diameter = C{11}     %Sink diameter
dist            = C{12}     %Src to Snk distance
tolerance       = C{13}     %Tolerance for distance

%%Load Source Data from id_dp_tm and Making it unique
dx = (sphere_diameter/2):sphere_diameter:(rlenx-sphere_diameter/2);
dy = (sphere_diameter/2):sphere_diameter:(rleny-sphere_diameter/2);
dz = (sphere_diameter/2):sphere_diameter:(rlenz-sphere_diameter/2);
xc_src = dx(nxi); yc_src = dy(nyi); zc_src = dz(nzi);

time_start = 0:(time_gap):(nfile-1-time_gap);
time_end = time_start + time_gap;
time_limit = nfile - 1 - time_end;

cond = time_limit > time_gap;
time_start = time_start(cond);
time_end = time_end(cond);
time_limit = time_limit(cond);

tic
disp('-----------------------------------')
disp('Load Source data')
[id_source,dp_source,tm_source,~,~,~] = load_source_data_v2(nxi,nyi,nzi,dpath);
toc
fprintf('\n')

%%Determining the locations of all Sinks at the chosen distance from Source
tic
disp('-----------------------------------')
disp('Locating Neighbours')
[sinks_loc] = neighbours(nxi, nyi, nzi, rlenx, rleny, rlenz, dist, sphere_diameter, tolerance);
%Sinks loc is a NX3 matrix with x,y,z location of all sinks at a distance
%dist from source
toc
fprintf('\n')

%%Cumulative variables
Nexit_cum       = cell(size(sinks_loc,1),numel(time_start));
Ntot_cum        = zeros(size(sinks_loc,1),numel(time_start));
Nhat_cum        = cell(size(sinks_loc,1),numel(time_start));
Nsfr_cum        = zeros(size(sinks_loc,1),(nfile - time_gap));
src_id_unique   = cell(numel(time_start),1);
src_dp_unique   = cell(numel(time_start),1);
src_tm_unique   = cell(numel(time_start),1);

for j=1:numel(time_start)
    tic
    disp('-----------------------------------')
    disp('Remove data from t = 0 to time_gap')
    [id_source0,dp_source0,tm_source0] = remove_data_till_tm(id_source,dp_source,...
        tm_source,time_start(j));
    toc
    fprintf('\n')

    tic
    disp('-----------------------------------')
    disp('Make data unique')
    [src_id_unique0,src_dp_unique0,src_tm_unique0] = make_unique(id_source0,dp_source0,tm_source0);
    toc
    fprintf('\n')

    %%Remove source particles after nfile-time_limit
    tic
    disp('-----------------------------------')
    disp('Removing data from t = time_limit to t = nfile -1')
    [src_id_unique{j},src_dp_unique{j},src_tm_unique{j}] = remove_data_after_tm(src_id_unique0,src_dp_unique0,src_tm_unique0,time_end(j));
    toc
end

%%Iterating over every Source Sink combination for given source
for i=1:size(sinks_loc,1)
    
    %%Load Sink data
    tic
    nx_snk = sinks_loc(i,1);
    ny_snk = sinks_loc(i,2);
    nz_snk = sinks_loc(i,3);
    text = strcat(['Sink: nxi ',num2str(nx_snk),...
        ', nyi ',num2str(ny_snk),...
        ', nzi ',num2str(nz_snk)]);
    disp(text);
    tic
    disp('-----------------------------------')
    disp('Load Sink data')

    [snk_id,snk_dp,snk_tm] = load_sink_data_v2(nx_snk,ny_snk,nz_snk,dpath);
    toc
    fprintf('\n')


    for j=1:numel(time_start)
        %%Identify particles within Source that enter Sink
        s = time_start(j)+1;
        e = time_end(j)+1;

        src_id = src_id_unique{j};
        src_tm = src_tm_unique{j};
        src_dp = src_dp_unique{j};
        tic
        disp('-----------------------------------')
        disp('Comparing time taken for particles to reach from Source to Sink')
        [src2snkDt] = src2snk_v02(src_id,src_tm,snk_id,snk_tm);
        [Nexit,Nhat,Ntot_cum(i,j)] = src2snkstats_v02(src_dp,src2snkDt,time_limit(j));

        Nhat = Nhat./Ntot_cum(i,j);
        tmp = zeros(size(Nhat));
        tmp((time_gap+1):end) = Nhat(1:(end - time_gap));
        Nsfr_cum(i,s:end) = Nsfr_cum(i,s:end) + Nhat - tmp;

        % plot(tmp); hold on;
        % plot(Nhat);
        % plot(Nsfr_cum(i,:));
        % % close all;
        toc
    end
end

%% Saving workspace
tic
disp('-----------------------------------')
disp('Saving the workspace')
text = strcat(['Source_nxi_',sprintf('%03d',round(nxi)),...
        '_nyi_',sprintf('%03d',round(nyi)),...
    '_nzi_',sprintf('%03d',round(nzi))]);
source_loc = [nxi, nyi, nzi];
save(fullfile(fpath, text),'Nsfr_cum','Ntot_cum','sinks_loc','time_start','time_end','time_limit')
disp('Workspace Saved')
toc
