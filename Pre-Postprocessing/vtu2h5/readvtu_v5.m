clc, close all, clear all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution file reader (readvtu)
infile = 'par00001.vtu';
sln_num = 6;
lrp_num = 12;
tag_num = 3;

[sln,lrp,tag,npart,time,step] = readvtu(infile,sln_num,lrp_num,tag_num);

%%
writevtu('parpuff00001.vtu',sln,lrp,tag,time,step);

%%
times = linspace(0,300,61);
for i=1:61
    aux = strcat('parpuff',pad(num2str(i,'%i'),5,'left','0'),'.vtu');
    files(i,:) = aux;
end
writepvd('par.pvd',files,times);
