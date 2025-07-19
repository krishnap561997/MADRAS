clc, clear all, close all
% V8: Process both par and parout without recompiling
%%
prompt = 'Number of Files?';
nfile = input(prompt)

%nfile = 20000;
%Variables - from Nek5000 simulation
ndiam = 10;
ninj = 2;	%ACH 2.5
nrank = 256;
nstep = 1250;	%ACH 5
%nstep = 625;	%ACH 2.5
%nstep = 2500;	%ACH 10
%max_tag03 = ninj*ndiam*4;
max_tag03 = round(1.5*ninj*ndiam*4); 

% Call readvtu
% sln_num = 6;
% lrp_num = 17;
% tag_num = 3;
% file_prefix = 'par';

 sln_num = 6;
 lrp_num = 18;
 tag_num = 5;
 file_prefix = 'parout';

for n=2:nfile
    infile = strcat(file_prefix,sprintf('%05d.vtu\n', n));
    tic
    [sln,lrp,tag,npart,time,step] = readvtu(infile,sln_num,lrp_num,tag_num);
    toc
    cond = tag(:,2)~=0;
    if sum(cond)==0
	continue
    end

    for i=1:size(sln,2)
        tmp = sln(:,i);
        sln2(:,i) = tmp(cond);
    end
    for i=1:size(lrp,2)
        tmp = lrp(:,i);
        lrp2(:,i) = tmp(cond);
    end
    for i=1:size(tag,2)
        tmp = tag(:,i);
        tag2(:,i) = tmp(cond);
    end
    
    sln = sln2;
    tag = tag2; 
    lrp = lrp2;
    tag_unique(:,1) = tag(:,2).*max_tag03 + tag(:,3);
    tag_time(:,1) = single((tag(:,2)+1)./nstep);

    % Save hdf5 files
    disp(['Save .h5 file for ',infile])
    tic
    outfile = strcat(file_prefix,sprintf('%05d.h5', n));
    h5create(outfile,'/g1/time',1,'Datatype','single')
    h5write(outfile,'/g1/time',time)
    h5create(outfile,'/g1/tagtime',size(single(lrp(:,2))),'Datatype','single')
    h5write(outfile,'/g1/tagtime',tag_time)
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

    tag2 = unique(tag_unique);
    counts = hist(tag_unique,tag2); counts = counts';
    rep2 = counts(counts~=1);
    rep = tag2(counts~=1);

    if numel(rep2)>0
       break;
    end

    
    % Clear arrays
     clear sln lrp tag npart time step tag_unique tag_time 
     clear sln2 lrp2 tag2 counts rep2 rep
end

clear sln lrp tag npart time step tag_unique sln_num lrp_num tag_num file_prefix outfile infile tag_time
clear sln2 lrp2 tag2 
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

for n=2:nfile
    infile = strcat(file_prefix,sprintf('%05d.vtu\n', n));
    tic
    [sln,lrp,tag,npart,time,step] = readvtu(infile,sln_num,lrp_num,tag_num);
    toc

%    if(n==1)
%        max_id1 = max(tag(:,1));
%        max_id3 = max(tag(:,3));
%    end
    cond = tag(:,2)~=0;
    for i=1:size(sln,2)
        tmp = sln(:,i);
        sln2(:,i) = tmp(cond);
    end
    for i=1:size(lrp,2)
        tmp = lrp(:,i);
        lrp2(:,i) = tmp(cond);
    end
    for i=1:size(tag,2)
        tmp = tag(:,i);
        tag2(:,i) = tmp(cond);
    end
    
    sln = sln2;
    tag = tag2; 
    lrp = lrp2;
    tag_unique(:,1) = tag(:,2).*max_tag03 + tag(:,3);
    tag_time(:,1) = single((tag(:,2)+1)./nstep);
    
    % Save hdf5 files
    disp(['Save .h5 file for ',infile])
    tic
    outfile = strcat(file_prefix,sprintf('%05d.h5', n));
    h5create(outfile,'/g1/time',1,'Datatype','single')
    h5write(outfile,'/g1/time',time)
    h5create(outfile,'/g1/tagtime',size(single(tag(:,2))),'Datatype','single')
    h5write(outfile,'/g1/tagtime',tag_time)
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
    tag2 = unique(tag_unique);
    counts = hist(tag_unique,tag2); counts = counts';
    rep2 = counts(counts~=1);
    rep = tag2(counts~=1);

    if numel(rep2)>0
       break;
    end

    
    % Clear arrays
     clear sln lrp tag npart time step tag_unique tag_time 
     clear sln2 lrp2 tag2 rep2 rep counts
end

