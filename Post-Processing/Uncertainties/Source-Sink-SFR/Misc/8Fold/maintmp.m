clc, clear all, close all;

%% Inputs
fileID = fopen('input.txt','r');
C = textscan(fileID,'%s %f %f %f %f','Delimiter','\n');

fpath   = char(strip(strip(C{1},'left',''''),'right',''''))
%File dump location. CHANGE IT!
spath   = pwd;     %id_dp_tm Location
fclose(fileID);

nxt             = C{2}
nyt             = C{3}
nzt             = C{4}
dist            = C{5}     %ADD IT TO INPUT.TXT GENERATING FUNCTION


for k=13:13  
    for i=k:13
        for j=3:3           
	    tic
            fname0 = strcat(['Source_nxi_',sprintf('%03d',round(i)),...
                '_nyi_',sprintf('%03d',round(j)),...
                '_nzi_',sprintf('%03d',round(k))]);
            disp(fname0);
            load(fullfile(fpath,fname0));
            nreal = size(Nexit_cum,2);

            nsnk = size(sinks_loc,1);
            src = [i,j,k];
            snk = sinks_loc;

            Nexit = cell(nsnk,nreal,8);
            Nhat = cell(nsnk,nreal,8);
            Ntot = zeros(nsnk,nreal,8);
            src_loc_c = zeros(nsnk,8);
            snk_loc_c = zeros(nsnk,8);


            for n=1:nsnk %Iterating over all sinks for corresponding 8 sources

                nx_locs = [src(1); snk(n,1)];                
                nz_locs = [src(3); snk(n,3)];
		
		if (src(1) == snk(n,1)) && (src(3) == snk(n,3))
                    [nx, nz, ~] = eightfoldsymmetry_locs(nx_locs(1), nz_locs(1), nxt, nzt, 1);
                    nx_src = nx(1,:);
                    ny_src = src(2);
                    nz_src = nz(1,:);
                    nx_snk = nx(1,:);
                    ny_snk = snk(n,2);
                    nz_snk = nz(1,:);
                else
                    [nx, nz, ~] = eightfoldsymmetry_locs(nx_locs, nz_locs, nxt, nzt, 2);
                    nx_src = nx(1,:);
                    ny_src = src(2);
                    nz_src = nz(1,:);
                    nx_snk = nx(2,:);
                    ny_snk = snk(n,2);
                    nz_snk = nz(2,:);
                end

                for m=1:8 %Iterating over 8 sources - 8 fold average
                    fname = strcat(['Source_nxi_',sprintf('%03d',round(nx_src(m))),...
                    '_nyi_',sprintf('%03d',round(ny_src)),...
                    '_nzi_',sprintf('%03d',round(nz_src(m)))]);
                    srctmp = encode_array([nx_src(m),ny_src,nz_src(m)],...
                            nxt,nyt,nzt,1);
                    load(fullfile(fpath,fname));
                    snktmp = encode_array([nx_snk(m),ny_snk,nz_snk(m)],...
                        nxt,nyt,nzt,1);
                    snkdata = encode_array(sinks_loc,...
                        nxt,nyt,nzt,1);     
                    [~, loc] = ismember(snktmp,snkdata);
                    % srctmp
                    % snktmp
                    src_loc_c(n,m) = srctmp;
                    snk_loc_c(n,m) = snktmp;
                    
                    for o=1:nreal
                        Nexit{n,o,m} = Nexit_cum{loc,o};
                        Nhat{n,o,m} = Nhat_cum{loc,o};
                        Ntot(n,o,m) = Ntot_cum(loc,o);
                        plot(Nhat{n,o,m}); hold on;
                    end
                    
                end
            end
            fname0 = strcat(['Source_8f_nxi_',sprintf('%03d',round(i)),...
                '_nyi_',sprintf('%03d',round(j)),...
                '_nzi_',sprintf('%03d',round(k)),'.mat']);
            save(fullfile(spath,fname0),'Ntot','Nhat','Nexit','src_loc_c','snk_loc_c','time_start',...
                'time_end','time_limit','time_end');
            disp('File Saved');
            toc
            clearvars -except fpath spath nxt nyt nzt dist i j k;
        end
    end    
end
