clear; clc; close all;

%% Add path
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %s %s %f %f %f %f','Delimiter','\n');

datapath    = char(strip(strip(C{1},'left',''''),'right',''''))
fpath       = char(strip(strip(C{2},'left',''''),'right',''''))
saveloc     = char(strip(strip(C{3},'left',''''),'right',''''))
%parent  = char(strip(strip(C{5},'left',''''),'right',''''));
fclose(fileID);
nx_src      = C{4}
ny_src      = C{5}
nz_src      = C{6}
ACH         = C{7}

fpath0 = strcat([datapath,'/Source_nxi_',...
                sprintf('%03d',round(nx_src)),'_nyi_',sprintf('%03d',round(ny_src)),...
                '_nzi_',sprintf('%03d',round(nz_src))]);
% nx_src = 3; ny_src = 5; nz_src = 7;
%Dd_cum = cell(25,8,25);

for nx=1:25
    for ny=1:8
        for nz=1:25
            tic
            fname = strcat(['Sink_nxi_',sprintf('%03d',round(nx)),...
                        '_nyi_',sprintf('%03d',round(ny)),...
                        '_nzi_',sprintf('%03d',round(nz)),'.mat']);
            disp(strcat(['nx: ',num2str(nx),'. ny: ',num2str(ny),'. nz: ',num2str(nz)]));
            load(fullfile(fpath0,fname));
            V_sink = (pi/6)*(0.4^3);
            n_src2snk0 = (Nhat./Ntot)./V_sink;
            t0 = 0:(size(Nhat,2)-1);
            clear Nexit Nhat Ntot 

            %%

            activity = 1;
            activity_list = {'Mouth_Breathing','Nose_Breathing','Singing','Soft_Singing',...
                'Speaking','Whispering'};

            rho = 1000; g = 9.81;   nu = 1.56e-5;   Height=3.2; % Height in BB is 12 ft=3.66m
            V_room = 10*10*3.2; mu = 1.81e-05;

            rc = sqrt(9*mu*ACH*Height/(7200*rho*g));    % Obtained from Bazant&Bush    
            %% Finding the parameters for ACH 5 - non dimensionalized settling velocity and concentration
            ACH0 = 5;
            unique_dp0 = [0.2 1 5 10 15 20 25 30 40 50];
            unique_rp0 = unique_dp0.*(1e-06)./2;

            unique_vs0   = 2*rho*(unique_rp0.^2)*g/(9*mu);   % With no Re correction
            Vs_tilde0 = unique_vs0./(Height.*ACH0/3600);

            %% Converting from ACH 5 to desired ACH

            unique_vs = (Height.*ACH0/3600).*Vs_tilde0;         %Converting vstilde to Vs for chosen ACH
            unique_rp = ((unique_vs.*9.*mu)./(2*rho)).^(0.5);    %Using Vs to convert to radius

            n_src2snk_tmp = n_src2snk0(1,:).*ones(size(n_src2snk0)).*ACH0./ACH;
            t_tmp = t0.*ACH0./ACH;
            t = 0:fix(max(t_tmp));

            for i=1:size(n_src2snk_tmp,1)
                n_src2snk1(i,:) = spline(t_tmp,n_src2snk_tmp(i,:),t);
            end

            [r, n_so,Qb] = choose_activity_HPG(activity,rc,fpath);

            for i=1:numel(t)
                n_src2snk = pchip(unique_rp, n_src2snk1(:,i), r);

                d_num_src2snk = n_src2snk.*n_so;
                d_den_src2snk = n_so;

                num_src2snk = trapz(r, d_num_src2snk);
                den_src2snk = trapz(r, d_den_src2snk);
                N_src2snk = num_src2snk./den_src2snk;

                nd(i) = max([N_src2snk,0]);
                if i==1
                    Dd(i) = nd(i);
                else
                    Dd(i) = Dd(i-1)+nd(i);
                end
            end
            %Dd_cum{nx,ny,nz} = Dd;
            Dd_cum(nx,ny,nz,:) = Dd;
            nd_cum(nx,ny,nz,:) = nd;
            clear 'num_src2snk' 'den_src2snk' 'N_src2snk' 'Dd' 'nd'
            toc
        end
    end
end
savename = strcat(['Source_nxi_',...
                sprintf('%03d',round(nx_src)),'_nyi_',sprintf('%03d',round(ny_src)),...
                '_nzi_',sprintf('%03d',round(nz_src)),'_ACH_',num2str(ACH),'.mat']);
save(fullfile(saveloc,savename),'Dd_cum','t','ACH','nd_cum');

a+b+c
