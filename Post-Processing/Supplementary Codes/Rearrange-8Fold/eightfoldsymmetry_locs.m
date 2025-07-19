 function [nx, nz, sym_locs] = eightfoldsymmetry_locs(nx_locs, nz_locs, nx_tot, nz_tot, n_locs)
% nx_locs, nz_locs are the x and z indices of the spheres.
% n_locs are the number of spheres we are considering. 
x = 1:nx_tot;
z = 1:nz_tot;
[X Z] = meshgrid(x,z);
locations = zeros(size(X));

nx = zeros(n_locs, 8);
nz = zeros(n_locs, 8);
% nx(:,1) = reshape(nx_locs, [n_locs, 1]);
% nz(:,1) = reshape(nz_locs, [n_locs, 1]);
for i=1:n_locs
    locations(nx_locs(i), nz_locs(i)) = i;
end

sym_locs = zeros(nx_tot,nz_tot,8);
sym_locs(:,:,1) = locations;
sym_locs(:,:,2) = fliplr(locations);
sym_locs(:,:,3) = flipud(locations);
sym_locs(:,:,4) = flipud(fliplr(locations));
sym_locs(:,:,5) = flipud(rot90(fliplr(locations),-1));
sym_locs(:,:,6) = rot90(rot90(fliplr(locations),-1),2); 
sym_locs(:,:,7) = rot90(fliplr(locations),-1); 
sym_locs(:,:,8) = fliplr(rot90(fliplr(locations),-1));

for i=1:8
    tmp = reshape(sym_locs(:,:,i),[nx_tot,nz_tot]);
    for j=1:n_locs
        cond = tmp == j;
        nx(j,i) = X(cond);
        nz(j,i) = Z(cond);
    end
end

% Figures to check
% X = (X-1).*0.4 + 0.2;
% Z = (Z-1).*0.4 + 0.2;
% for i=1:8
%     subplot(2,4,i)
%     contour(X,Z,sym_locs(:,:,i),n_locs);
%     xlim([0 10]); ylim([0 10]);
%     set(gca,'FontSize',24)
%     set(gca,'linewidth',3)
% end
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 16, 8], 'PaperUnits', 'Inches', 'PaperSize', [16, 8]);