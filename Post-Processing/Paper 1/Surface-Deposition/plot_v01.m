function [] = plot_v01(x_s,y_s,z_s,unique_diam,unique_surf,x,y,z,fpath,fontsize)

%% Plot number density of deposit on floor for largest particle
% Now construct the histogram

[xf,zf] = meshgrid(x,z);
[xw,yw] = meshgrid(x,y);

for i=1:numel(unique_diam)
    for j=1:numel(unique_surf)
        xz{i,j} = [x_s{i,j}, z_s{i,j}];
        xy{i,j} = [y_s{i,j}, x_s{i,j}];
        yz{i,j} = [y_s{i,j}, z_s{i,j}];
    end
end

for i=1:numel(unique_diam)
    for j=1:3
        if(j==1) % Floor
            counts = hist3(xz{i,j}, size(xf));
            counts_avg = (1/8)*(counts+fliplr(counts)...
             + flipud(counts+fliplr(counts))...
             + flipud(rot90(fliplr(counts),-1))+rot90(rot90(fliplr(counts),-1),2)...
             + rot90(fliplr(counts),-1)+fliplr(rot90(fliplr(counts),-1))...
               );            
            figure((i-1)*numel(unique_diam)+j)
            figure_2d_v01(counts_avg,xf,zf,unique_diam(i),'x','z','Floor',fpath,fontsize)

        elseif(j==3) % Ceiling
            counts = hist3(xz{i,j}, size(xf));
            counts_avg = (1/8)*(counts+fliplr(counts)...
             + flipud(counts+fliplr(counts))...
             + flipud(rot90(fliplr(counts),-1))+rot90(rot90(fliplr(counts),-1),2)...
             + rot90(fliplr(counts),-1)+fliplr(rot90(fliplr(counts),-1))...
               );            
            figure((i-1)*numel(unique_diam)+j)
            figure_2d_v01(counts_avg,xf,zf,unique_diam(i),'x','z','Ceiling',fpath,fontsize)
            
        else % sidewalls
            counts = hist3(yz{i,2}, size(xw)) + hist3(yz{i,4}, size(xw))...
                +    hist3(xy{i,5}, size(xw)) + hist3(xy{i,6}, size(xw));
            counts_avg = (1/8)*(counts+fliplr(counts));            
            figure((i-1)*numel(unique_diam)+j)
            figure_2d_v01(counts_avg,xw,yw,unique_diam(i),'x,z','y','Wall',fpath,fontsize)
        end
    end
end
