clear all; close all; clc;
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %f %f %f %f','Delimiter','\n');
spath           = char(strip(strip(C{1},'left',''''),'right',''''))       %id_dp_tm location

K               = C{2}
ACH             = C{3}
dist            = C{4}
tol		= C{5}
ctr=0;

% K = 0.04548;
% ACH = 5;
dx = 0.2:0.4:9.8;
dz = dx;
dy = 0.2:0.4:3;
L = 10;
H = 3.2;
r = round(0.001:0.001:20,3);
t = 0:1600;
[Y,X,Z] = meshgrid(dy,dx,dz);
[T, R] = meshgrid(t,r);
Rm = exp(-(R.^2)./(4.*K.*T));
% tol = 0.02;
ctr = 0;
R_sum = zeros(size(t));
C0 = (1./(4.*pi.*K.*t).^1.5).*exp(-ACH.*t./3600);
for i=1:numel(dx)
    for j=1:numel(dy)
        for k=1:numel(dz)
            x0 = dx(i);
            y0 = dy(j);
            z0 = dz(k);
            disp([i,j,k]);
            d = round(sqrt((X-x0).^2 + (Y-y0).^2 + (Z-z0).^2),3);
%            cond = (d <= dist+tol) | (d >= dist-tol);
            cond = round(d,1) == round(dist,1);
            [~, loc] = ismember(dist, r);
            R0 = Rm(loc, :);
            x_snk = X(cond);
            y_snk = Y(cond);
            z_snk = Z(cond);
            for l = 1:numel(x_snk)
                x = x_snk(l);
                y = y_snk(l);
                z = z_snk(l);
%                 scatter3(x0,z0,y0,'r');
%                 hold on;
%                 scatter3(x,z,y,'b');
%                 scatter3(-x0,z,y,'b');
%                 xlim([0 10]);
%                 ylim([0 10]);
%                 zlim([0 3.2]);
                d_x1 = min(round(sqrt((x+x0).^2 + (y-y0).^2 + (z-z0).^2),3),max(r));
                [~, locx1] = ismember(d_x1, r);
                d_y1 = min(round(sqrt((x-x0).^2 + (y+y0).^2 + (z-z0).^2),3),max(r));
                [~, locy1] = ismember(d_y1, r);
                d_z1 = min(round(sqrt((x-x0).^2 + (y-y0).^2 + (z+z0).^2),3),max(r));
                [~, locz1] = ismember(d_z1, r);
                d_x2 = min(round(sqrt((x+x0-2*L).^2 + (y-y0).^2 + (z-z0).^2),3),max(r));
                [~, locx2] = ismember(d_x2, r);
                d_y2 = min(round(sqrt((x-x0).^2 + (y+y0-2*H).^2 + (z-z0).^2),3),max(r));
                [~, locy2] = ismember(d_y2, r);
                d_z2 = min(round(sqrt((x-x0).^2 + (y-y0).^2 + (z+z0-2*L).^2),3),max(r));
                [~, locz2] = ismember(d_z2, r);
                
                tmp = R0 + Rm(locx1,:) + Rm(locy1,:) + Rm(locz1,:) + ...
                    Rm(locx2,:) + Rm(locy2,:) + Rm(locz2,:);
                R_sum = R_sum + tmp;
                ctr = ctr+1;
                Nhat_EDM(ctr,:)=tmp;
            end   
        end
    end
end

Nhat_mean = C0.*R_sum./ctr;
Nhat_EDM = C0.*Nhat_EDM;

text = strcat('EDMData_Dist_',num2str(dist),'m.mat');
save(fullfile(spath,text),'R_sum','C0','ctr','dist','K','ACH',...
    'Nhat_mean','Nhat_EDM','-v7.3');
disp('Data Saved');
% fname= strcat(['EDM_dist_',num2str(round(dist,1)),'m.mat']);
% save(fname,'n_EDM','n_EDM2','C0','R0','R_sum');
