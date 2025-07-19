function [x_s,y_s,z_s,unique_diam,unique_surf,x,y,z] = sort_v01(dp,time,surf,loc)

dx = 0.1; % Must be consistent with round(loc,!)

%% Sort particles by diameter and surface of deposit
unique_diam = unique(dp);
unique_surf = unique(surf);

%% Define the meshgrid on each of the 6 surfaces
x = 0:dx:max(loc(:,1));
y = 0:dx:max(loc(:,2));
z = 0:dx:max(loc(:,3));

loc_rd = round(loc,1); % Rounding must be consistent with dx

for i=1:numel(unique_diam)
    for j=1:numel(unique_surf)
        x_s{i,j}   = loc_rd(dp==unique_diam(i) & surf==unique_surf(j),1);
        y_s{i,j}   = loc_rd(dp==unique_diam(i) & surf==unique_surf(j),2);
        z_s{i,j}   = loc_rd(dp==unique_diam(i) & surf==unique_surf(j),3);
    end
end

