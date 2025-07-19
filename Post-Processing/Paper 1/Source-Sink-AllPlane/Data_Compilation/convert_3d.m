function [count_spheres] = convert_3d(Nhat_end, sinks_loc)
rlenx = 10;
rleny = 3.2;
rlenz = 10;
d = 0.4;
dx = d/2:d:(rlenx - d/2);
dy = d/2:d:(rleny - d/2);
dz = d/2:d:(rlenz - d/2);
nx = numel(dx);
ny = numel(dy);
nz = numel(dz);

data = Nhat_end;

[midY, midX, midZ] = meshgrid(dy, dx, dz);

count_spheres = zeros(size(midX,1),size(midY,2),size(midZ,3),size(data,1));

for i=1:size(sinks_loc,1)
    xsnk = sinks_loc(i,:);
    count_spheres(xsnk(1), xsnk(2),xsnk(3),:) = data(:,i);
end
