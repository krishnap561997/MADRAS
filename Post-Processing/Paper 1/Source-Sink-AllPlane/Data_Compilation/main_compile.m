nxi = 13;
nyi = 3;
nzi = 13;
rlenx = 10;
rleny = 3.2;
rlenz = 10;
sphere_diameter = 0.4;
nx = 25;
ny = 8;
nz = 25;

X = find_symm(nxi, nzi, rlenx, rlenz, sphere_diameter);

pFolder = '/blue/bala1s/krishnap.kalivel/COVID/Sou2Snk_AllPlane/room_30';
txt1 = strcat(['Source_nxi_',sprintf('%03d',round(nxi)),...
    '_nyi_',sprintf('%03d',round(nyi)),...
    '_nzi_',sprintf('%03d',round(nzi))]);

pFolder = fullfile(pFolder, txt1)


for i=1:size(X,1)
    txt2 = strcat(['Source_nxi_',sprintf('%03d',round(X(i,1))),...
    '_nyi_',sprintf('%03d',round(nyi)),...
    '_nzi_',sprintf('%03d',round(X(i,2)))]);
    pFolder2 = fullfile(pFolder,txt2)
    Nhat_end_c = [];
    sinks_loc_c = [];
    
    for i1=1:nx
	for j=1:ny
		for k=1:nz
		        fname = strcat(['Sink_nxi_',sprintf('%03d',round(i1)),...
		        '_nyi_',sprintf('%03d',round(j)),...
		        '_nzi_',sprintf('%03d',round(k)),'.mat'])
			fname2 = fullfile(pFolder2,fname);
		        load(fname2);
		        processdata;
        		sinks_loc = [i1,j,k];
		        Nhat_end_c = [Nhat_end_c'; Nhat_end']';
		        sinks_loc_c = [sinks_loc_c; sinks_loc];
		end
	end
    end
    
    Nhat_end = Nhat_end_c;
    sinks_loc = sinks_loc_c;
    count_spheres = convert_3d(Nhat_end_c, sinks_loc_c);
    save(fullfile(pFolder, txt2), 'count_spheres','Nhat_end','sinks_loc')
end
                            
                           
