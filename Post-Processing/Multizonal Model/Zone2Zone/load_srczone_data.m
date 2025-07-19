function [id_src, dp_src, loc_src] = load_srczone_data(fpath, start_t, Z_src, rlenx, rlenz, outlet_c)

n = start_t + 1;
tic
infile = strcat(sprintf('par%05d.h5\n', n));
disp('-----------------------------------')
disp(['Reading ',infile])
time0  = h5read(fullfile(fpath,infile),'/g1/time');% single
disp(strcat(['t = ',num2str(time0),' s']))
id0    = h5read(fullfile(fpath,infile),'/g1/tag'); % unit32
dp0    = h5read(fullfile(fpath,infile),'/g1/dp');  % single
loc0   = h5read(fullfile(fpath,infile),'/g1/loc'); % single


[id_src, dp_src, loc_src] = load_zone_data(Z_src, id0, dp0, loc0, rlenx, rlenz, outlet_c);