function count = get_zonal_stats(fpath, start_t, i, rlenx, rlenz, outlet_c, id_src, unique_dp)

count = zeros(numel(unique_dp), 4);

n = start_t + i;
tic
infile = strcat(sprintf('par%05d.h5\n', n));
disp('-----------------------------------')
disp(['Reading ',infile])
time0  = h5read(fullfile(fpath,infile),'/g1/time');% single
disp(strcat(['t = ',num2str(time0),' s']))
id0    = h5read(fullfile(fpath,infile),'/g1/tag'); % unit32
dp0    = h5read(fullfile(fpath,infile),'/g1/dp');  % single
loc0   = h5read(fullfile(fpath,infile),'/g1/loc'); % single
toc

for j=1:4
    tic
    disp('-----------------------------------')
    disp(['Zone ',num2str(j)])
    [id_snk0, dp_snk0, ~] = load_zone_data(j, id0, dp0, loc0, rlenx, rlenz, outlet_c);
    for k=1:numel(unique_dp)        
        id_snk = id_snk0(dp_snk0 == unique_dp(k));
        count(k,j) = sum(ismember(id_src{k},id_snk));
    end
    toc
end

a=1;