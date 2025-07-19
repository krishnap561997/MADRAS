function [A] = find_symm(nxi, nzi, rlenx, rlenz, sphere_diameter)

x1 = 0; x2 = 10;
z1 = 0; z2 = 10;
d = sphere_diameter;

dx = d/2:d:(rlenx-d/2);
dz = d/2:d:(rlenz-d/2);
x = dx(nxi); z = dz(nzi);

A = [nxi, nzi];
wlx = round(min(abs(x-x1), abs(x-x2)),2);
wlz = round(min(abs(z-z1), abs(z-z2)),2);

for i=1:numel(dx)
    for k=1:numel(dz)
        x0 = dx(i);
        z0 = dz(k);
        wlx0 = round(min(abs(x0-x1), abs(x0-x2)),2);
        wlz0 = round(min(abs(z0-z1), abs(z0-z2)),2);
        cond1 = (wlx==wlx0) && (wlz==wlz0);
        cond2 = (wlx==wlz0) && (wlz==wlx0);
        if cond1 || cond2
            tmp = [i,k];
            A = [A; tmp];
        end
    end
end

A = unique(A,'rows');

