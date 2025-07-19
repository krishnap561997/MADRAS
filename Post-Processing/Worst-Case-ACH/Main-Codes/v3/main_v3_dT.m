clear; clc; close all;
%% dT1 = dT2
%% Add path
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %f %f %f %f %f %f %f','Delimiter','\n');

fpath   = char(strip(strip(C{1},'left',''''),'right',''''))
fclose(fileID);
d       = C{2}
T0     = C{3}
t_min  = C{4}
t_max  = C{5}
ACH_eff_min = C{6}
ACH_eff_max = C{7}
eta_ef  = C{8} 
lam_v   = 0.3/3600

activity_list = {'Mouth_Breathing','Nose_Breathing','Singing','Soft_Singing',...
    'Speaking','Whispering'};
activity = 1;
H = 3.2;
V = 320;

ACH_min = min(ACH_eff_min,round(ACH_eff_min./eta_ef))
ACH_max = max(ACH_eff_max,round(ACH_eff_max./eta_ef))


dT1_list = (t_min:0.05:t_max).*60;
ACH_list = ACH_min:0.1:ACH_max;
T0 = T0*60;


for i=1:numel(ACH_list)
    ACH = ACH_list(i);
    ACS = ACH/3600;
    tic
    % Get well-mixed concentration
    n_wm = get_nwm(ACH,H,V,eta_ef,lam_v,activity,2.*max(dT1_list)+T0);
    % Get distance specific concentration
    n_d = get_nd(n_wm,ACH,eta_ef,d);

    for j=1:numel(dT1_list)
        % Cumulative Exposure at Sink
        dT1 = dT1_list(j);
        disp(strcat(['ACH: ',num2str(ACH),'. dT1: ',num2str(T0/60),' min']));
        D_d(i,j) = get_exposure(n_d,ACS,V,dT1,T0,dT1);
    end
    toc
end
dT1_list = dT1_list./60;
T0 = T0/60;
dT = 2.*dT1_list;

fname = strcat(['CumExp_dist_',num2str(d),...
    'm_T0_',num2str(T0),'_eff_',num2str(eta_ef),'_dTequal.mat']);
save(fullfile(fpath,fname),'D_d','d','dT1_list','dT','T0','ACH_list','V','H');

a+b+c
