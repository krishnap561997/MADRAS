clear all; close all; clc;
%location to store data, id_dp_tm location, file size, tlim, room
%dimensions, sphere diameter, distance between source and sink
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %f %f %f %f %f %f','Delimiter','\n');

fpath   = char(strip(strip(C{1},'left',''''),'right',''''))        %Save location

dist            = C{2}      %Distance
mc_nreal        = C{3}      %Number of Trials
time_gap        = C{4}      %dt
i               = C{5}
j               = C{6}
k               = C{7}

%% Codes
%% Load Data
% fpath = fullfile(pwd,strcat(['dist_',num2str(dist),'m']));

fname1 = strcat(['Source_r31_nxi_',sprintf('%03d',round(i)),...
                '_nyi_',sprintf('%03d',round(j)),...
                '_nzi_',sprintf('%03d',round(k))]);
fname2 = strcat(['Source_r34_nxi_',sprintf('%03d',round(i)),...
                '_nyi_',sprintf('%03d',round(j)),...
                '_nzi_',sprintf('%03d',round(k))]);

load(fullfile(fpath,fname1));
ctr = 1;
for loc=1:size(snk_loc_c,1)
    output = encode_array(snk_loc_c(loc,:),25, 8, 25, 2);
    if output(1,2) == j        
        ctr = 0;
        break;
    end
end

if ctr == 0
    loc = 1;
end

nreal = size(Nexit_cum,2);
nsnk = size(Nexit_cum,3);
tlim = size(Nexit_cum,4);
DFR_conc = reshape(Nexit_cum(loc,:,:,:)./Ntot_cum(loc,:,:),[nreal.*nsnk,tlim]);
SFR_conc = reshape(Nhat_cum(loc,:,:,:)./Ntot_cum(loc,:,:),[nreal.*nsnk,tlim]);


load(fullfile(fpath,fname2));
tmp = reshape(Nexit_cum(loc,:,:,1:tlim)./Ntot_cum(loc,:,:),[nreal.*nsnk,tlim]);
DFR_conc = [DFR_conc; tmp];
tmp = reshape(Nhat_cum(loc,:,:,1:tlim)./Ntot_cum(loc,:,:),[nreal.*nsnk,tlim]);
SFR_conc = [SFR_conc; tmp];

clear Nexit_cum Nhat_cum Ntot_cum i j k nreal nsnk loc output fname1 fname2 ...
    time_end time_limit time_start src_loc_c snk_loc_c tmp;

%% SFR
close all;
index_start = 1:time_gap:(tlim);
index_end = tlim - index_start + 1;
tmp0 = zeros(size(SFR_conc));
tmp0(:,(time_gap+1):tlim) = -SFR_conc(:,1:(tlim-time_gap));
SFR_tgap = SFR_conc + tmp0;
SFR_ind = zeros(size(SFR_conc,1),size(SFR_conc,1),numel(index_start));
for i=1:numel(index_start)
    SFR_ind(:,index_start(i):tlim,i) = SFR_tgap(:,1:index_end(i));
end
t = 0:(tlim-1);
% figure(1)
% subplot(1,2,1)
% for i=1:size(SFR_conc,1)
%     plot(SFR_tgap(i,:),'LineWidth',2); hold on;
% end
% xlabel('$t$','interpreter','latex');
% ylabel('$\int_0 ^t c_\delta(t) dt$ - $\int_0^{t - \Delta t} c_\delta (t - \Delta t) dt$',...
%     'Interpreter','latex')
% box on;
% xlim([0 1610]);
% set(gca,'FontSize',15);
% set(gca,'LineWidth',2);
% 
% subplot(1,2,2)
% for i=1:numel(index_start)
%     plot(SFR_ind(1,:,i),'LineWidth',2); hold on;
% end
% box on;
% xlim([0 1610]);
% set(gca,'FontSize',15);
% set(gca,'LineWidth',2);
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 18, 6], 'PaperUnits', 'Inches', 'PaperSize', [18 6]);
% 
% name = strcat(['Figure0_trials_',num2str(mc_nreal),'_dt_',num2str(time_gap),'.png']);
% print('-dpng','-r300',fullfile(fpath,name));

%% Randomize
close all;
mc_rno = randi([1 size(DFR_conc,1)], mc_nreal, numel(index_start));
SFR_approx = zeros(mc_nreal,tlim);
tic
for i=1:mc_nreal
    disp(strcat(['Monte-Carlo Trial Number: ',num2str(i)]));
    for j=1:numel(index_start)        
        SFR_approx(i,:) = SFR_approx(i,:) + SFR_ind(mc_rno(i,j),:,j);
    end
end
toc
%% Plots
close all;
figure(1)
nlines = 10;
index = randi([1, mc_nreal], nlines, 1);

subplot(1,2,1)
for i=1:nlines
    plot(t,SFR_approx(index(i),:),'LineWidth',2); hold on;
end

xlabel('$t$','interpreter','latex');
ylabel('SFR approximation',...
    'Interpreter','latex')
box on;
xlim([0 1610]);
set(gca,'FontSize',15);
set(gca,'LineWidth',2);

subplot(1,2,2)
tmp = SFR_approx(:,end);
avg = round(mean(tmp),3);
sigma = round(std(tmp,0),3);
histogram(SFR_approx(:,end),'Normalization','pdf');
xlabel('SFR approximation ($t = 1610 s$)','Interpreter','latex');
ylabel('PDF');
xline(avg,'LineWidth',2);
set(gca,'FontSize',15);
set(gca,'LineWidth',2);
text = strcat([num2str(mc_nreal),' trials. ',...
    '$\langle C \rangle$ = ',num2str(avg)...
    ,'. (var$(C))^{1/2}$ = ',num2str(sigma)]);
sgtitle(text,'interpreter','latex','FontSize',20);

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 18, 6], 'PaperUnits', 'Inches', 'PaperSize', [18 6]);

name = strcat(['MC_real_',num2str(mc_nreal),'.png']);
print('-dpng','-r300',fullfile(fpath,name))

fname = strcat(['SFRData_dist_',num2str(dist),'m_trials_',num2str(mc_nreal),'_dt_',num2str(time_gap),'s.mat']);
save(fullfile(fpath,fname),'SFR_approx');

a+b+c
