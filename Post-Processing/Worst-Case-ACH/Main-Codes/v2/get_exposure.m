function D_d = get_exposure(n_d,ACS,V_room,T_si,T_so,dT)

n_d_t = V_room*ACS*n_d;
t =  0:(numel(n_d)-1);
t_t = t.*ACS;

if dT<=0
    T_si = T_si + dT;
else
    T_so = T_so - dT;
end

if T_si <= T_so
    cond = t<=T_si;
    tmp1 = n_d_t(cond);
    tmp2 = t_t(cond);
    plot(t(cond),tmp1./(V_room*ACS^2)); hold on;
    D_d = trapz(tmp2,tmp1)./(V_room*ACS^2);
elseif T_si > T_so
    cond = t<=T_si & t>=(T_si-T_so);
    tmp1 = n_d_t(cond);
    tmp2 = t_t(cond);
    plot(t(cond),tmp1./(V_room*ACS^2)); hold on;
    D_d = trapz(tmp2,tmp1)./(V_room*ACS^2);
end
