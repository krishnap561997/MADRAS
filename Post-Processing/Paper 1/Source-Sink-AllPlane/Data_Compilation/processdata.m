Nhat_c = [];
Ntot_c = [];
Nhat_c = [];
Nhat_c_end = [];
Nhat_end = [];
d = 0.4;
V_sink = (1/6)*pi*(d^3);

%Average Concentration for ACH 10
%lambda2 = [1.2197    1.2185    1.1803    0.9825    0.7429    0.5499    0.4081    0.3086    0.1870    0.1205];

%Average concentration for ACH 2.5
%lambda2 = [4.8194    4.8180    3.9107    2.2206    1.2554    0.7732    0.5089    0.3477    0.1738    0.0981];

%Average Concentration for ACH 5
%lambda2 = [2.4380    2.4378    2.2287    1.5792    1.0389    0.6955    0.4845    0.3527    0.1968    0.1130];

%Area under curve
lambda2 = [2.17774603928958 2.12367132678614 1.83688021574628 1.34081255359887 0.923173751287389 0.640483123652415 0.458480546220332 0.339666553718651 0.203138443849159 0.127402296997310]
unique_dp = [0.2 1 5 10 15 20 25 30 40 50];


Nhat_c = Nhat./Ntot.*ones(size(Nhat));
Nhat_c_end = Nhat_c(:,end);
Nhat_end = Nhat_c_end./(V_sink.*lambda2'.*ones(size(Nhat_c_end)));

