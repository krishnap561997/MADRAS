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

function [xsphere,ysphere,zsphere] = center_location(nxi,nzi,sphere_diameter,rlenx,rleny,rlenz)
% nxi is an integer between 1 and rlenx/sphere_diameter, which in this case is 1 to 25
% nzi is an integer between 1 and rlenz/sphere_diameter, which in this case is 1 to 25

nyt = round(rleny / sphere_diameter);

% rlenx = 10.0;   % Domain size in x
% rleny =  3.2;   % Domain size in y
% rlenz = 10.0;   % Domain size in z

%sphere_diameter = 0.4; % Must be chosen so that:
% rlenx/nx = rlenz/nz = rleny = ny = sphere_diameter. nx, ny, and nz are the number of spheres in x, y, and z.
x = sphere_diameter/2 : sphere_diameter : rlenx - sphere_diameter/2;
z = sphere_diameter/2 : sphere_diameter : rlenz - sphere_diameter/2;
%y = sphere_diameter/2 : sphere_diameter : rleny - sphere_diameter/2;

ysphere = (sphere_diameter/2 : sphere_diameter : rleny - sphere_diameter/2)'; 
xsphere = x(nxi) * ones(size(ysphere));
zsphere = z(nzi) * ones(size(ysphere));

disp('-----------------------------------')
disp('Source center locations: x, y, z')
[xsphere(:) ysphere(:) zsphere(:)] 
fprintf('\n')

