fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %s %f %f %f','Delimiter','\n');
par     = char(strip(strip(C{1},'left',''''),'right',''''))
parout  = char(strip(strip(C{2},'left',''''),'right',''''))
fclose(fileID);

addpath(genpath(parout))
addpath(genpath(par))
nfile	= C{3}
t_endinj= C{4}
ndp     = C{5};


id_unq = cell(ndp,1);
Ntot = zeros(ndp,1);
for n=2:(nfile)
    infile = strcat(sprintf('parout%05d.h5\n', n));
    disp('-----------------------------------')
    disp(['Reading ',infile])
    id_1    = h5read(infile,'/g1/tag'); % unit32
    dp_1    = h5read(infile,'/g1/dp');  % single
    if numel(unique(dp_1))==ndp
        unique_dp = unique(dp_1);
    end
    for i=1:ndp
        [tf idx] = ismember(dp_1, unique_dp(i));
        tmp = id_1(tf);
        if numel(tmp)==0
           continue;
        end
        if n==2
            tmp2 = tmp;
        else
            tmp2 = [id_unq{i}; tmp];
        end
        id_unq{i} = tmp2;
    end
end

for i=1:10
    tmp = unique(id_unq{i});
    text = strcat(['Non Unique: dp: ',num2str(unique_dp(i)),' microns. N: ',num2str(numel(id_unq{i}))]);
    disp(text);
    text = strcat(['Unique: dp: ',num2str(unique_dp(i)),' microns. N: ',num2str(numel(tmp))]);
    disp(text);
    Ntot(i,1) = numel(id_unq{i});
end

