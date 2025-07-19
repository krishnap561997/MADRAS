clear; clc; close all;
%% Add path
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %f','Delimiter','\n');

fpath   = char(strip(strip(C{1},'left',''''),'right',''''))
fclose(fileID);
d     = C{2}

eta_ef = 1;                 % Filter efficiency
lam_v = 0.3/3600;           % Viral deactivation rate (From Bazant & Bush)
activity_list = {'Mouth_Breathing','Nose_Breathing','Singing','Soft_Singing',...
    'Speaking','Whispering'};
activity = 1;
H = 3.2; V = 320;           % Height and Volume of the room (m, m^3)
dT1_list = 0:2:40;
%dT1_list = 0:0.5:2;
dT2_list = 0:2:40;
%dT2_list = 0:0.5:2;
T0_list = 1:1:30;
%T0_list = 1:0.5:5;
%% Low ACH
ACH = 0.5;

for i=1:numel(T0_list)
    for j=1:numel(dT1_list)
        for k=1:numel(dT2_list)
            tic
            dT1 = dT1_list(j)*60;
            dT2 = dT2_list(k)*60;
            T0 = T0_list(i)*60;
            T_si = T0*60 + dT2;
            T_so = T0*60 + dT1;

            disp(strcat(['Low ACH. dt1: ',num2str(dT1/60),'. T0: ',num2str(T0/60),...
                '. dT2: ',num2str(dT2/60)]));
            ACS = ACH/3600;
                       
            n_wm = get_nwm(ACH, H,V, eta_ef, lam_v, activity, max(T_si+dT1,T_so));
            %% Concentration for specified source-sink separation
            n_d = get_nd(n_wm, ACH, eta_ef,d,dT1);
  
            %% Cumulative Exposure at Sink
            D_d_low(i,j,k) = get_exposure(n_d,ACS,V,T_si,T_so,dT1);
            toc
        end
    end
end


%% High ACH
ACH = 10;

for i=1:numel(T0_list)
    for j=1:numel(dT1_list)
        for k=1:numel(dT2_list)
            dT1 = dT1_list(j)*60;
            dT2 = dT2_list(k)*60;
            T0 = T0_list(i)*60;
            T_si = T0*60 + dT2;
            T_so = T0*60 + dT1;
            tic
            disp(strcat(['High ACH. dt1: ',num2str(dT1/60),'. T0: ',num2str(T0/60),...
                '. dT2: ',num2str(dT2/60)]));
            ACS = ACH/3600;
                       
            n_wm = get_nwm(ACH, H,V, eta_ef, lam_v, activity, max(T_si+dT1,T_so));
            %% Concentration for specified source-sink separation
            n_d = get_nd(n_wm, ACH, eta_ef,d,dT1);
 
            %% Cumulative Exposure at Sink
            D_d_high = get_exposure(n_d,ACS,V,T_si,T_so,dT1);
            toc
        end
    end
end

fname = strcat(['HighvsLow_','d_',num2str(d),'.mat']);
save(fullfile(fpath,fname),'D_d_high','D_d_low','T0_list','dT1_list','dT2_list');

a+b+c
