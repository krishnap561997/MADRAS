% data_vel and below are not present in the old parout format
% room_24_nolan has the old parout format. 
% Therefore, comment from data_vel and below and recompile
% if rerunning the room_24_nolan case
function [id,dp,time,surf,loc] = deposit_v01(nfile)

for n=2:nfile
    tic
    infile = strcat(sprintf('parout%05d.h5\n', n));
    disp('-----------------------------------')
    disp(['Reading ',infile])
    outfile     = strcat(sprintf('parout%05d.h5', n));
    data_time   = h5read(outfile,'/g1/time');
    data_id     = h5read(outfile,'/g1/tag');
    data_dp     = h5read(outfile,'/g1/dp');
    data_loc    = h5read(outfile,'/g1/loc');
    data_cycle  = h5read(outfile,'/g1/cyc');
    data_wall   = h5read(outfile,'/g1/wll');
    data_vel   = h5read(outfile,'/g1/vel'); % Particle velocity from the solution array (only useful for dp>10mu, it is 0 otherwise)
    data_uel   = h5read(outfile,'/g1/uel'); % Fluid velocity at marker
    data_uelp  = h5read(outfile,'/g1/uelp'); % Langevin perturbation at marker
    data_tot   = h5read(outfile,'/g1/utot'); % Fluid + Langevin at marker
    data_lan   = h5read(outfile,'/g1/lan'); % Some langevin related data

    toc
    fprintf('\n')
% -------------------------------------------------------------------------
%% Extract info of deposited particles
    tic
    if n==2
        id      = data_id;
        dp      = data_dp;
        loc     = data_loc;
        surf    = data_wall;
        time    = data_dp*0 + data_time;
    else
        id   = [id    ; data_id];
        dp   = [dp    ; data_dp];
        loc  = [loc   ; data_loc];
        surf = [surf  ; data_wall];  
        time = [time  ; data_dp*0 + data_time];        
    end
    
    disp([num2str(numel(data_id)),' particles have exited the domain'])
    toc
    fprintf('\n')
% -------------------------------------------------------------------------
end

