function [id_1,dp_1,time,x,y,z,exit_surf] = record_times_v02(nfile)

% clc, clear all
n = 1;
tic
infile = strcat(sprintf('par%05d.h5\n', n));
disp('-----------------------------------')
disp(['Reading ',infile])
time_1  = h5read(infile,'/g1/time');% single
id_1    = h5read(infile,'/g1/tag'); % unit32
dp_1    = h5read(infile,'/g1/dp');  % single
loc_1   = h5read(infile,'/g1/loc'); % single
toc
% -------------------------------------------------------------------------
 
x = loc_1(:,1);
y = loc_1(:,2);
z = loc_1(:,3);

exit_surf = 0*single(dp_1);

time = 0*single(id_1);
% -------------------------------------------------------------------------
for n=2:nfile
%     for n=4
    tic
    infile = strcat(sprintf('par%05d.h5\n', n));
    disp('-----------------------------------')
    disp(['Reading ',infile])
    time_n  = h5read(infile,'/g1/time');% single
    id_n    = h5read(infile,'/g1/tag'); % unit32
    dp_n    = h5read(infile,'/g1/dp');  % single
    loc_n   = h5read(infile,'/g1/loc'); % single
    toc
% -------------------------------------------------------------------------
%% Record the exit time of particles for each output
    tic
    disp('-----------------------------------')
    disp('Finding particles that have exited the domain')
    [Lia,~] = ismember(id_1,id_n);
    time(Lia==0 & time==0) = time_n;    % Lia==0  -> Particles that now exited
                                        % time==0 -> Particles previously in domain
    toc
    
    tic
    disp('-----------------------------------')
    disp('Recording exit location of particles')
    % Particles deposit on floor
    exit_surf(time==time_n & y<0.1) = 1; 
    
    % Particles exit through outlet
    exit_surf(time==time_n & y>3 & abs(x-5)<0.3 & abs(z-5)<0.3) = 2;    
    toc
% -------------------------------------------------------------------------
    
%% Update the x,y,z location of particles still ariborne
    tic
    disp('-----------------------------------')
    disp('Updating x,y,z location of airborne particles')
    [tf,idx] = ismember(id_n,id_1); % Find particles that are still ariborne
    x(idx(idx~=0)) = loc_n(tf,1);   % Update x location of airborne particles
    y(idx(idx~=0)) = loc_n(tf,2);   % Update y location of airborne particles
    z(idx(idx~=0)) = loc_n(tf,3);   % Update z location of airborne particles
        
    disp([num2str(sum(Lia==0)),' particles exited the domain'])
    toc
end
% end
