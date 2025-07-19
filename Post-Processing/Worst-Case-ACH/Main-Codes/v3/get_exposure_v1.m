function D_d = get_exposure(n_d,ACS,V_room,dT1,T0,dT2)

n_d((dT1+T0+dT2+2):end)=[];

t =  -dT1:(T0+dT2);
n_d(t<0)=[];
t(t<0)=[];
n_d_t = V_room*ACS*n_d;
t_t = t.*ACS;

cond = t<=(T0+dT2) & t >= dT2;
tmp1 = n_d_t(cond);
tmp2 = t_t(cond);
D_d = trapz(tmp2,tmp1)./(V_room*ACS^2);

