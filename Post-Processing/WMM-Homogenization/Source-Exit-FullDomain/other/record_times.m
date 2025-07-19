function [data_id_1,data_dp_1,data_time] = record_times(nfile)

% clc, clear all
n = 1;
tic
infile = strcat(sprintf('par%05d.h5\n', n));
disp('-----------------------------------')
disp(['Reading ',infile])
time_1        = h5read(infile,'/g1/time');% single
data_id_1     = h5read(infile,'/g1/tag'); % unit32
data_dp_1     = h5read(infile,'/g1/dp');  % single
data_loc_1    = h5read(infile,'/g1/loc'); % single

toc
% -------------------------------------------------------------------------
  
data_time = 0*single(data_id_1);
for n=2:nfile
   tic
   infile = strcat(sprintf('par%05d.h5\n', n));
   disp('-----------------------------------')
   disp(['Reading ',infile])
   time        = h5read(infile,'/g1/time');% single
   data_id     = h5read(infile,'/g1/tag'); % unit32
   data_dp     = h5read(infile,'/g1/dp');  % single
   data_loc    = h5read(infile,'/g1/loc'); % single

   disp('Finding particles that have exited the domain')
   [Lia,~] = ismember(data_id_1,data_id);
       
   data_time(Lia==0 & data_time==0) = time;

   disp([num2str(sum(Lia==0)),' particles exited the domain'])
   toc
end
