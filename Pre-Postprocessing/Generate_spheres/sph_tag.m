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

function [sph_id_cum,sph_dp_cum,sph_tm_cum] = sph_tag(sphere_diameter,xsphere,ysphere,zsphere,nfile)

%% Source geometric properties
num_ctr = numel(ysphere); % Number of spheres considered

%% Read data for all particles in the domain
for n=1:nfile
    tic
    infile = strcat(sprintf('par%05d.h5\n', n));
    disp('-----------------------------------')
    disp(['Reading ',infile])
    time        = h5read(infile,'/g1/time');% single
    data_id     = h5read(infile,'/g1/tag'); % unit32
    data_dp     = h5read(infile,'/g1/dp');  % single
    data_loc    = h5read(infile,'/g1/loc'); % single
    toc
% -------------------------------------------------------------------------
    xs = data_loc(:,1);
    ys = data_loc(:,2);
    zs = data_loc(:,3);
    %% Find particles inside ALL num_ctr spheres
    tic
    sph_cond    = cell(num_ctr,1);
    sph_id      = cell(num_ctr,1);
    sph_diam    = cell(num_ctr,1);
    sph_time    = cell(num_ctr,1);
    for l=1:num_ctr
        sph_cond{l} = (xs-xsphere(1)).^2 ...
                    + (ys-ysphere(l)).^2 ...
                    + (zs-zsphere(l)).^2 <= (sphere_diameter/2)^2;
        sph_id{l}   = data_id(sph_cond{l});
        sph_diam{l} = data_dp(sph_cond{l});
        sph_time{l} = time*ones(size(sph_diam{l}));
    end
    disp('-----------------------------------')
    disp(['Finding particles inside the ',num2str(num_ctr),' spheres'])
    toc
    tic
    if(n==1)
        for l=1:num_ctr
            sph_id_cum{l} = sph_id{l};
            sph_dp_cum{l} = sph_diam{l};
            sph_tm_cum{l} = sph_time{l};
        end
    else
        for l=1:num_ctr
            sph_id_cum{l} = [sph_id_cum{l} ; sph_id{l}];
            sph_dp_cum{l} = [sph_dp_cum{l} ; sph_diam{l}];
            sph_tm_cum{l} = [sph_tm_cum{l} ; sph_time{l}];
        end
    end
    disp('-----------------------------------')
    disp(['Concatenating the ',num2str(num_ctr),' source matrices'])
    toc
end
% -------------------------------------------------------------------------

end

