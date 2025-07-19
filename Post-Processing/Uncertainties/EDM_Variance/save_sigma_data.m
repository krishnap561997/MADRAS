clear; clc; close all;
d=  [0.4,1,2,6,8];
fpath = '/blue/bala/krishnap.kalivel/COVID';
for i=1:numel(d)
    tic
    disp(d(i));
    fname = strcat(['EDMData_Dist_',num2str(d(i)),'m.mat']);
    load(fullfile(fpath,fname));
    Nhat_EDM(:,1) = 0; Nhat_mean(1) = 0;
    N_EDM = zeros(size(Nhat_EDM));
    N_mean = zeros(size(Nhat_mean));
    for j=1:size(Nhat_EDM,2)
        N_EDM(:,j) = trapz(Nhat_EDM(:,1:j),2);
        N_mean(j) = trapz(Nhat_mean(:,1:j));
        j
    end
    sigma2_EDM = std(N_EDM,1,1);
    mean_EDM = N_mean;
    sname = strcat(['StatsEDM_dist_',num2str(d(i)),'m.mat']);
    save(sname,'ACH','K','sigma2_EDM','mean_EDM','dist');
    toc
end
