function [Nwall,Nfloor,Noutlet,Ntot,unique_dp,unique_time] = exit_time_v02(dp,time,exit_surf)

%% Compute the number of particle sizes involved
unique_dp   = unique(dp);
unique_time = unique(time);
delta_out   = min(diff(unique_time));
text1 = ['Time difference b/w outputs is ',num2str(delta_out),' sec'];
text2 = ['Final output time that is read is ',num2str(max(unique_time)),' sec'];
disp(text1)
disp(text2)

% figure
% frame_h = get(handle(gcf),'JavaFrame');
% set(frame_h,'Maximized',1);

Ntot    = zeros(numel(unique_dp),1);
Nwall   = zeros(numel(unique_dp),numel(unique_time));
Nfloor  = zeros(numel(unique_dp),numel(unique_time));
Noutlet = zeros(numel(unique_dp),numel(unique_time));

for i=1:numel(unique_dp)
    for j = 1:numel(unique_time)

        idx         = (dp==unique_dp(i));
        Ntot(i,1)   = sum(idx);        

        idx         = (dp==unique_dp(i) & (time>0) & (time<=unique_time(j))...
                        & exit_surf==0);
        Nwall(i,j)  = sum(idx);       

        idx         = (dp==unique_dp(i) & (time>0) & (time<=unique_time(j))...
                        & exit_surf==1);
        Nfloor(i,j)  = sum(idx);      

        idx         = (dp==unique_dp(i) & (time>0) & (time<=unique_time(j))...
                        & exit_surf==2);
        Noutlet(i,j)  = sum(idx);  
    end
end
    