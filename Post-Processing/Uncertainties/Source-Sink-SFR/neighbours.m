%% *** MAtlab turbulent Dispersion Room-scale Analysis Software (MADRAS)*** 
% 
% This software is for calculating turbulent dispersal statistics, such as
% concentration and flux, from Lagrangian particle tracks generated from
% DNS, LES, and RANS.
% 
% The software has also been specialized for predicting pathogen concentrat
% -ion in indoor spaces to compute quantities such as cumulative exposure
% time (CET) and Safe Occupancy Limit
% 
% The use/distribution of this code is allowed ONLY after approval by Prof.
% S. Balachandar *(bala1s@ufl.edu)*, the project PI. This code was written
% and conceptualized at the University of Florida at the Computational
% Multiphysics Group and the Center for Compressible Multiphase Turbulence
% by Prof. S. Balachandar, Prof. Nadim Zgheib, Dr. Jorge Salinas, and K.A.
% Krishnaprasad.
% *************************************************************************

function [sinks_loc] = neighbours(nxi, nyi, nzi, rlenx, rleny, rlenz, dist, sphere_diameter, tolerance)
%Spheres Center locations - x, y, z
x_c = sphere_diameter/2:sphere_diameter:(rlenx - sphere_diameter/2);
y_c = sphere_diameter/2:sphere_diameter:(rleny - sphere_diameter/2);
z_c = sphere_diameter/2:sphere_diameter:(rlenz - sphere_diameter/2);

sinks_loc = [];
source_loc = [x_c(nxi), y_c(nyi), z_c(nzi)];

nx = numel(x_c);
ny = numel(y_c);
nz = numel(z_c);

l = 1;

for i=1:nx
    for j=1:ny
        for k =1:nz
            temp = [x_c(i), y_c(j), z_c(k)];
            distance = round(norm(source_loc - temp),2);
            if distance <= dist + tolerance && distance >= dist - tolerance
                loc = [i, j, k];
                sinks_loc = [sinks_loc; loc];
            end
        end
    end
end




