d = 8;
fpath = strcat(['/blue/bala/krishnap.kalivel/COVID/HigherACHBetter/dist',...
    num2str(d),'m/Data_dT1']);
dT2_list = 0:45;
T0_list = 1:45;
for i1 = 1:45
    for i2 = 0:45
        T0 = i1;
        dT2 = i2;
        fname = strcat(['CumExp_dist_',num2str(d),...
            'm_T0_',num2str(T0),'min_dT2_',num2str(dT2),'min.mat'])
        load(fullfile(fpath,fname));
        ACH_listnew = min(ACH_list):0.01:max(ACH_list);
        dT1_listnew = min(dT1_list):0.01:max(dT1_list);
        D_d(isnan(D_d)) = 0;
        [X,Y] = meshgrid(dT1_list,ACH_list);
        Z = D_d;

        [Xn,Yn]=meshgrid(dT1_listnew,ACH_listnew);
        Zn = interp2(X,Y,Z,Xn,Yn,'makima');
        for i=1:numel(dT1_listnew)
            exposure = Zn(:,i);
            [max_val, loc] = max(exposure);
            ACHw(i,i1,i2+1) = ACH_listnew(loc);
            half_val = max_val/2;
            quart_val = max_val/4;
            [~, loc2] = min(abs(exposure(1:loc)-half_val));
            ACHw_half(i,1,i1,i2+1) = ACH_listnew(loc2);
            [~, loc3] = min(abs(exposure(loc:end)-half_val));
            ACHw_half(i,2,i1,i2+1) = ACH_listnew(loc3+loc-1);


            [~, loc4] = min(abs(exposure(1:loc)-quart_val));
            ACHw_quart(i,1,i1+1,i2+1) = ACH_listnew(loc4);
            [~, loc5] = min(abs(exposure(loc:end)-quart_val));
            ACHw_quart(i,2,i1+1,i2+1) = ACH_listnew(loc5+loc-1);
        end

    end
end
fname = strcat(['CumExp_dT1_dist_',num2str(d),'m.mat']);
save(fname, 'ACHw_quart','ACHw_half','ACHw','dT1_listnew','ACH_listnew',...
    'T0_list','dT2_list','d');

a+b+c

