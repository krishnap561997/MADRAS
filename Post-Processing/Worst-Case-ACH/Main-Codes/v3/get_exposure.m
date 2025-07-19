function D_d = get_exposure(n_d,ACS,V_room,dT1,T0,dT2)

n_d = n_d(1:round(dT1+T0+dT2+1));
n_d_t = V_room*ACS*n_d;
t =  0:round(dT1+T0+dT2);
t_t = t.*ACS;

tmp1 = n_d_t;
tmp2 = t_t;
D_d_all = trapz(tmp2,tmp1)./(V_room*ACS^2);
D_d_dT1 = 0;
D_d_dT2 = 0;
if dT1 ~=0
    cond = t<=dT1;
    tmp1 = n_d_t(cond);
    tmp2 = t_t(cond);
    D_d_dT1 = trapz(tmp2,tmp1)./(V_room*ACS^2);
end

if dT2~=0
    cond = t<=dT2;
    tmp1 = n_d_t(cond);
    tmp2 = t_t(cond);
    D_d_dT2 = trapz(tmp2,tmp1)./(V_room*ACS^2);
end
D_d = D_d_all - D_d_dT1 - D_d_dT2;



