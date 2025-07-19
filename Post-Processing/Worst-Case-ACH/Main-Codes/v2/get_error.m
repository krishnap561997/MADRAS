function error = get_error(D_d,ACH,dist,T0)

modeldata = D_d;

fname = strcat(['T0_10min_ACH_',num2str(ACH),'_dist_',num2str(dist),'m.mat']);
load(fname);
rawdata = D_d(T0+1);
error = 100.*(modeldata-rawdata)./rawdata;