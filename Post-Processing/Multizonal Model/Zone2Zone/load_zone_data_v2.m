function [id, dp, loc] = load_zone_data_v2(zone, id0, dp0, loc0, rlenx, rleny, rlenz, outlet_c, dw)

rlx = 0; rly =0; rlz = 0; rrx = rlenx; rry = rleny; rrz =rlenz; 

cond_wall = (loc0(:,1) > rlx+dw) & (loc0(:,1) < rrx - dw) & (loc0(:,2) > rly + dw) & (loc0(:,2) < rry - dw) ...
	& (loc0(:,3) > rlz + dw) & (loc0(:,3) < rrz - dw);

id0 = id0(cond_wall);
dp0 = dp0(cond_wall);
loc0 = loc0(cond_wall,:);

disp('Removing particles close to wall')
disp(strcat(['X: Min: ',num2str(min(loc0(:,1))),'. Max:',num2str(max(loc0(:,1)))]));
disp(strcat(['Y: Min: ',num2str(min(loc0(:,2))),'. Max:',num2str(max(loc0(:,2)))]));
disp(strcat(['Z: Min: ',num2str(min(loc0(:,3))),'. Max:',num2str(max(loc0(:,3)))]));

switch zone
    case 1
        cond = (loc0(:,1) < outlet_c(1)) & ((loc0(:,3)-rlz)./(loc0(:,1) - rlx) > ((outlet_c(2)-rlz)./(outlet_c(1)-rlx))) ...
            & ((loc0(:,3)-rrz)./(loc0(:,1) - rlx) < ((outlet_c(2)-rrz)./(outlet_c(1)-rlx)));
    case 2
        cond = (loc0(:,3) < outlet_c(2)) & ((loc0(:,3)-rlz)./(loc0(:,1) - rlx) < ((outlet_c(2)-rlz)./(outlet_c(1)-rlx))) ...
            & ((loc0(:,3)-rlz)./(loc0(:,1) - rrx) > ((outlet_c(2)-rlz)./(outlet_c(1)-rrx)));
    case 3
        cond = (loc0(:,1) > outlet_c(1)) & ((loc0(:,3)-rrz)./(loc0(:,1) - rrx) > ((outlet_c(2)-rrz)./(outlet_c(1)-rrx))) ...
            & ((loc0(:,3)-rlz)./(loc0(:,1) - rrx) < ((outlet_c(2)-rlz)./(outlet_c(1)-rrx)));
    case 4
        cond = (loc0(:,3) > outlet_c(2)) & ((loc0(:,3)-rrz)./(loc0(:,1) - rrx) < ((outlet_c(2)-rrz)./(outlet_c(1)-rrx))) ...
            & ((loc0(:,3)-rrz)./(loc0(:,1) - rlx) > ((outlet_c(2)-rrz)./(outlet_c(1)-rlx)));
end

id = id0(cond);
dp = dp0(cond);
loc = loc0(cond,:);
