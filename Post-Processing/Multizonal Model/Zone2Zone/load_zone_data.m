function [id, dp, loc] = load_zone_data(zone, id0, dp0, loc0, rlenx, rlenz, outlet_c)

rlx = 0; rlz = 0; rrx = rlenx; rrz =rlenz;

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