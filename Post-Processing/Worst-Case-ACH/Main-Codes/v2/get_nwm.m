function n_wm = get_nwm(ACH, H,V, eta_ef, lam_v ,activity,T)

rho = 1000; g = 9.81; mu = 1.81e-05;            % Constants

rc          = sqrt(9*mu*ACH*H/(7200*rho*g));    % Obtained from Bazant&Bush
unique_rp = 0:(rc/10):rc;
[r, n_so,~] = choose_activity(activity,rc);
unique_vs   = 2*rho*(r.^2)*g/(9*mu);

lambda_wm = eta_ef*ACH/3600 + unique_vs/H  + lam_v;
for i=0:T
    n_room = (1-exp(-lambda_wm*i))./(lambda_wm*V);
    d_num_room = n_room.*n_so;
    d_den_room = n_so;
    num_room = trapz(r, d_num_room);
    den_room = trapz(r, d_den_room);
    n_wm(i+1) = num_room./den_room;
end