clc, clear all, close all
% V8: Process both par and parout without recompiling
%%
prompt = 'Number of Files?';
nfile = input(prompt)

%nfile = 20000;
%Variables - from Nek5000 simulation
npart = 7000;
ndiam = 10;
nrank = 256;
%nstep = 1250;	%ACH 5
nstep = 625;	%ACH 2.5
%nstep = 2500;	%ACH 10
max_id3 = npart*ndiam;
max_tag0 = round(npart*ndiam*nrank,-6);
buffer = 1e06; 

% Call readvtu
% sln_num = 6;
% lrp_num = 17;
% tag_num = 3;
% file_prefix = 'par';

 sln_num = 6;
 lrp_num = 18;
 tag_num = 5;
 file_prefix = 'parout';

for n=1:nfile
    infile = strcat(file_prefix,sprintf('%05d.vtu\n', n));
    tic
    [sln,lrp,tag,npart,time,step] = readvtu(infile,sln_num,lrp_num,tag_num);
    toc

%    if(n==1)
%        max_id1 = max(tag(:,1));
%        max_id3 = max(tag(:,3));
%    end
%    tag_unique(:,1) = tag(:,1).*max_id3 + tag(:,3);
    tag_unique(:,1) = (tag(:,1).*max_id3) ...
                + (tag(:,2)~=0).*((1+round(tag(:,2)./nstep)).*buffer + max_tag0)...
                + tag(:,3); 
    % Save hdf5 files
    disp(['Save .h5 file for ',infile])
    tic
    outfile = strcat(file_prefix,sprintf('%05d.h5', n));
    h5create(outfile,'/g1/time',1,'Datatype','single')
    h5write(outfile,'/g1/time',time)
    h5create(outfile,'/g1/tag',size(uint32(tag_unique(:,1))),'Datatype','uint32')
    h5write(outfile,'/g1/tag',uint32(tag_unique(:,1)))
    h5create(outfile,'/g1/dp',size(single(lrp(:,2))),'Datatype','single')
    h5write(outfile,'/g1/dp',single(lrp(:,2)))     
    h5create(outfile,'/g1/loc',size(single(sln(:,1:3))),'Datatype','single')
    h5write(outfile,'/g1/loc',single(sln(:,1:3)))
    
    if contains(file_prefix, 'parout')
        h5create(outfile,'/g1/cyc',size(uint32(tag(:,4))),'Datatype','uint32')
        h5write(outfile,'/g1/cyc',uint32(tag(:,4)))
        h5create(outfile,'/g1/wll',size(uint32(tag(:,5))),'Datatype','uint32')
        h5write(outfile,'/g1/wll',uint32(tag(:,5)))    
        %
        h5create(outfile,'/g1/vel',size(single(sln(:,4:6))),'Datatype','single')
        h5write(outfile,'/g1/vel',single(sln(:,4:6)))
        h5create(outfile,'/g1/uel',size(single(lrp(:,4:6))),'Datatype','single')
        h5write(outfile,'/g1/uel',single(lrp(:,4:6))) 
        h5create(outfile,'/g1/uelp',size(single(lrp(:,7:9))),'Datatype','single')
        h5write(outfile,'/g1/uelp',single(lrp(:,7:9)))
        h5create(outfile,'/g1/utot',size(single(lrp(:,10:12))),'Datatype','single')
        h5write(outfile,'/g1/utot',single(lrp(:,10:12)))
        h5create(outfile,'/g1/lan',size(single(lrp(:,13:16))),'Datatype','single')
        h5write(outfile,'/g1/lan',single(lrp(:,13:16)))    
    end
    
    toc
    
    % Clear arrays
     clear sln lrp tag npart time step tag_unique 
end

clear sln lrp tag npart time step tag_unique sln_num lrp_num tag_num file_prefix outfile infile
%clear all close all clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nfile = 20000;

% Call readvtu
 sln_num = 6;
 lrp_num = 17;
 tag_num = 3;
 file_prefix = 'par';

% sln_num = 6;
% lrp_num = 16;
% tag_num = 5;
% file_prefix = 'parout';

for n=1:nfile
    infile = strcat(file_prefix,sprintf('%05d.vtu\n', n));
    tic
    [sln,lrp,tag,npart,time,step] = readvtu(infile,sln_num,lrp_num,tag_num);
    toc

%    if(n==1)
%        max_id1 = max(tag(:,1));
%        max_id3 = max(tag(:,3));
%    end

    tag_unique(:,1) = (tag(:,1).*max_id3) ...
                + (tag(:,2)~=0).*((1+round(tag(:,2)./nstep)).*buffer + max_tag0)...
                + tag(:,3); 
    % Save hdf5 files
    disp(['Save .h5 file for ',infile])
    tic
    outfile = strcat(file_prefix,sprintf('%05d.h5', n));
    h5create(outfile,'/g1/time',1,'Datatype','single')
    h5write(outfile,'/g1/time',time)
    h5create(outfile,'/g1/tag',size(uint32(tag_unique(:,1))),'Datatype','uint32')
    h5write(outfile,'/g1/tag',uint32(tag_unique(:,1)))
    h5create(outfile,'/g1/dp',size(single(lrp(:,2))),'Datatype','single')
    h5write(outfile,'/g1/dp',single(lrp(:,2)))     
    h5create(outfile,'/g1/loc',size(single(sln(:,1:3))),'Datatype','single')
    h5write(outfile,'/g1/loc',single(sln(:,1:3)))
    
    if contains(file_prefix, 'parout')
        h5create(outfile,'/g1/cyc',size(uint32(tag(:,4))),'Datatype','uint32')
        h5write(outfile,'/g1/cyc',uint32(tag(:,4)))
        h5create(outfile,'/g1/wll',size(uint32(tag(:,5))),'Datatype','uint32')
        h5write(outfile,'/g1/wll',uint32(tag(:,5)))    
        %
        h5create(outfile,'/g1/vel',size(single(sln(:,4:6))),'Datatype','single')
        h5write(outfile,'/g1/vel',single(sln(:,4:6)))
        h5create(outfile,'/g1/uel',size(single(lrp(:,4:6))),'Datatype','single')
        h5write(outfile,'/g1/uel',single(lrp(:,4:6))) 
        h5create(outfile,'/g1/uelp',size(single(lrp(:,7:9))),'Datatype','single')
        h5write(outfile,'/g1/uelp',single(lrp(:,7:9)))
        h5create(outfile,'/g1/utot',size(single(lrp(:,10:12))),'Datatype','single')
        h5write(outfile,'/g1/utot',single(lrp(:,10:12)))
        h5create(outfile,'/g1/lan',size(single(lrp(:,13:16))),'Datatype','single')
        h5write(outfile,'/g1/lan',single(lrp(:,13:16)))    
    end
    
    toc
    
    % Clear arrays
     clear sln lrp tag npart time step tag_unique 
end
