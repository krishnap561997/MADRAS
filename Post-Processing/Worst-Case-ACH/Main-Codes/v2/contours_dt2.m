clear; clc; close all;
%% Add path
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %f','Delimiter','\n');

fpath   = char(strip(strip(C{1},'left',''''),'right',''''))
fclose(fileID);
T0     = C{2}

eta_ef = 1;                 % Filter efficiency
lam_v = 0.3/3600;           % Viral deactivation rate (From Bazant & Bush)
activity_list = {'Mouth_Breathing','Nose_Breathing','Singing','Soft_Singing',...
    'Speaking','Whispering'};
activity = 1;
H = 3.2; V = 320;           % Height and Volume of the room (m, m^3)
ACH_list = 0.3:0.1:10;
dist_list = 4:0.5:10;
dT2_list = 0:1:30;

for i=1:numel(ACH_list)
    ACH = ACH_list(i);
    ACS = ACH/3600;
    for j=1:numel(dist_list)
        d = dist_list(j);
        disp(strcat(['ACH :',num2str(ACH),'. Dist: ',num2str(d)]));
        tic
        for k=1:numel(dT2_list)
            dT = 0;
            T_si = T0*60 + dT2_list(k)*60;
            T_so = T0*60;
                       
            n_wm = get_nwm(ACH, H,V, eta_ef, lam_v, activity, max(T_si+dT,T_so));
            %% Concentration for specified source-sink separation
            n_d = get_nd(n_wm, ACH, eta_ef,d,dT);
    
            %% Cumulative Exposure at Sink
            D_d(i,j,k) = get_exposure(n_d,ACS,V,T_si,T_so,dT);
        end
        toc
    end
end

[X,Y] = meshgrid(dist_list,ACH_list);

fname = strcat(['ContourData_T0_',num2str(T0),'_',activity_list{activity},'dT2.mat']);
save(fullfile(fpath,fname),'X','Y','D_d','ACH_list','dist_list','T0','dT2_list');
