clear; clc; close all;
%% Add path
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %f %f %f %f %f %f %f %f','Delimiter','\n');

fpath   = char(strip(strip(C{1},'left',''''),'right',''''))
fclose(fileID);
T0      = C{2}
dT1     = C{3}
dT2     = C{4}
d_min   = C{5}
d_max   = C{6}
ACH_eff_min = C{7}
ACH_eff_max = C{8}
eta_ef  = C{9}


%% Data
lam_v   = 0.3/3600;
ACH_min = ACH_eff_min/eta_ef;
ACH_max = ACH_eff_max/eta_ef;
activity_list = {'Mouth_Breathing','Nose_Breathing','Singing','Soft_Singing',...
    'Speaking','Whispering'};
activity = 1;
H = 3.2;
V = 320;

T0 = T0.*60;
ACH_list = ACH_min:0.1:ACH_max;
dT1 = dT1*60;
dT2 = dT2*60;
d_list = d_min:0.1:d_max;


for i=1:numel(ACH_list)   
    ACH = ACH_list(i);
    ACS = ACH/3600;
    tic
    % Get well-mixed concentration
    n_wm = get_nwm(ACH,H,V,eta_ef,lam_v,activity,dT1+T0+dT2); 
      
    
    for j=1:numel(d_list)
        % Get distance specific concentration
        d = d_list(j);
        n_d = get_nd(n_wm,ACH,eta_ef,d);
        % Cumulative Exposure at Sink
        %T0 = T0_list(j);
        disp(strcat(['ACH: ',num2str(ACH),'. T0: ',num2str(T0/60),' min']));
        D_d(i,j) = get_exposure(n_d,ACS,V,dT1,T0,dT2);
    end
    toc
end
T0 = T0/60;
dT1 = dT1/60;
dT2 = dT2/60;

fname = strcat(['CumExp_T0_',num2str(T0),'min_eta_',num2str(eta_ef),...
    '_dT1_',num2str(dT1),'min_dT2_',num2str(dT2),'min.mat']);
save(fullfile(fpath,fname),'D_d','d_list','dT1','dT2','T0','ACH_list','V','H');

a+b+c
