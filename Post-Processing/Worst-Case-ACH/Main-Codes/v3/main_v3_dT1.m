clear; clc; close all;
%% Add path
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %f %f %f %f %f %f %f','Delimiter','\n');

fpath   = char(strip(strip(C{1},'left',''''),'right',''''))
fclose(fileID);
d       = C{2}
T0     = C{3}
dT2     = C{4}
t_min  = C{5}
t_max  = C{6}
ACH_min = C{7}
ACH_max = C{8}

%% Data
eta_ef  = 1
lam_v   = 0.3/3600

activity_list = {'Mouth_Breathing','Nose_Breathing','Singing','Soft_Singing',...
    'Speaking','Whispering'};
activity = 1;
H = 3.2;
V = 320;

dT1_list = (t_min:0.1:t_max).*60;
ACH_list = ACH_min:0.1:ACH_max;
T0 = T0*60;
dT2 = dT2*60;


for i=1:numel(ACH_list)
    ACH = ACH_list(i);
    ACS = ACH/3600;
    tic
    % Get well-mixed concentration
    n_wm = get_nwm(ACH,H,V,eta_ef,lam_v,activity,max(dT1_list)+T0+dT2);
    % Get distance specific concentration
    n_d = get_nd(n_wm,ACH,eta_ef,d);

    for j=1:numel(dT1_list)
        % Cumulative Exposure at Sink
        dT1 = dT1_list(j);
        disp(strcat(['ACH: ',num2str(ACH),'. dT1: ',num2str(T0),' min']));
        D_d(i,j) = get_exposure(n_d,ACS,V,dT1,T0,dT2);
    end
    toc
end
dT1_list = dT1_list./60;
T0 = T0/60;
dT2 = dT2/60;

fname = strcat(['CumExp_dist_',num2str(d),...
    'm_T0_',num2str(T0),'min_dT2_',num2str(dT2),'min.mat']);
save(fullfile(fpath,fname),'D_d','d','dT1_list','dT2','T0','ACH_list','V','H');

a+b+c
