function n_d = get_nd(n_wm, ACH, eta_ef,d,dT)

dist_list = [0.4,1,2,4,6,8,10];
lwr = [0, 0.33, 0.66];
upr = [0.33, 0.66, 1];
lwr_d = dist_list(1:end-1);
upr_d = dist_list(2:end);
cond = eta_ef > lwr & eta_ef <= upr;
cond_d = d > lwr_d & d <= upr_d;

fpath = '/home/krishnap.kalivel/MATLAB_src/Higher-ACH-Better/Main-Codes/v2/Data';
name = strcat([fpath,'/','Gamma_ACH_5_eff_',num2str(lwr(cond)),'.mat']);
load(name);
NaNchk = isnan(gamma_extrap);
gamma_extrap(NaNchk)=0;
gamma_lwr = gamma_extrap;
name = strcat([fpath,'/','Gamma_ACH_5_eff_',num2str(upr(cond)),'.mat']);
load(name);
NaNchk = isnan(gamma_extrap);
gamma_extrap(NaNchk)=0;
gamma_upr = gamma_extrap;


wt = (eta_ef-lwr(cond))/(upr(cond)-lwr(cond));
gamma1 = wt.*gamma_upr + (1-wt).*gamma_lwr;

wt_d = (d-lwr_d(cond_d))/(upr_d(cond_d)-lwr_d(cond_d));
pos = find(cond_d==1);
gamma_dlwr = gamma1(pos,:);
gamma_dupr = gamma1(pos+1,:);

gamma2 = round(wt_d,3).*gamma_dupr + round(1-wt_d,3).*gamma_dlwr;

t0 = 0:(numel(gamma2)-1);
t_adj = t0.*5/ACH;
t = 0:round(max(t_adj));
gamma = interp1(t_adj,gamma2,t);

if numel(n_wm) <= numel(gamma)
    gamma = gamma(1:numel(n_wm));
else
    gamma = [gamma, gamma(end).*ones(1,numel(n_wm)-numel(gamma))];
end
n_d = gamma.*n_wm;


% Added to v2. Cuts off the first dt seconds of data from n_d if source
% appears before the sink
if dT>0
    n_d(1:(dT))=[];
end
    
