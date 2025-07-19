function [Nexit,Ntot,unique_dp,unique_time] = exit_time(dp,time)

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
Nexit   = zeros(numel(unique_dp),numel(unique_time));

for i=1:numel(unique_dp)
    for j = 1:numel(unique_time)

        idx         = (dp==unique_dp(i));
        Ntot(i,1)   = sum(idx);        

        idx         = (dp==unique_dp(i) & (time>0) & (time<=unique_time(j)));
        Nexit(i,j)  = sum(idx);  
    end
end
    
%         subplot(3,1,1)
% %         plot(toto(1:end-20),100*Nexit(i,1:end-20)./Ntot(i,1:end-20),...
%         semilogx(unique_time,100*Nexit(i,:)./Ntot(i,:),...
%             '-','LineWidth',5)
%         hold on
%         ylabel('% particles that exit')
%         set(gca,'FontSize',15)
% % 
%         subplot(3,1,2)
%         semilogy(unique_time,Ntot(i,:),'-','LineWidth',5)
%         ylabel('N_{tot}')
%         set(gca,'FontSize',15)
%         hold on
% 
%         subplot(3,1,3)
%         semilogy(unique_time,Nexit(i,:),'-','LineWidth',5)
%         hold on
% 
%         Legend{i} = strcat(['d_p =',' ',num2str(unique_dp(i)*1e6),' \mum']);   
%         legend(Legend);
%         xlabel('exit time (s)')
%         ylabel('N_{exit}')
%         set(gca,'FontSize',15)
% end
% exit_mat = 1e9*ones(numel(distribution(:,1)),numel(unique_time)-1);
% 
% A = repmat(distribution(:,1),[1 length(unique_dp)]);
% [~,closestIndex] = min(abs(A-unique_dp'));
% closestValue = distribution(closestIndex,1);
% disp(['particles are of size ',num2str(1e6*closestValue'),' microns'])
% 
% exit_mat(closestIndex,:) = Nexit./Ntot;
% exit_mat(exit_mat==1e9)=nan;
% 
% for i1=1:numel(exit_mat(1,:))
%     nanx = isnan(exit_mat(:,i1));
%     t    = 1:numel(distribution(:,1));
%     exit_mat(nanx,i1) = interp1(t(~nanx), exit_mat(~nanx,i1), t(nanx));
% end
% 
% %% Total ejected cough volume
% Np_inj = repmat(distribution(:,2),1,numel(exit_mat(1,:)));
% dp_inj = repmat(distribution(:,1),1,numel(exit_mat(1,:)));
% Vp = sum(Np_inj.*(pi/6).*dp_inj.^3.*(1-exit_mat),1);
% Np = sum(Np_inj.*(1-exit_mat),1);
% 
% %% Plot of residence time
% figure
% frame_h = get(handle(gcf),'JavaFrame');
% set(frame_h,'Maximized',1);
% subplot(2,1,1)
% plot(unique_time(1:end-1),Vp/Vp(1),'o','MarkerSize',5,'LineWidth',3)
% % xlabel('time since ejection (s)')
% ylabel('Normalized particle volume')
% set(gca,'FontSize',15)
% subplot(2,1,2)
% plot(unique_time(1:end-1),Np/Np(1),'x','MarkerSize',5,'LineWidth',3)
% xlabel('time since ejection (s)')
% ylabel('Normalized particle count')
% set(gca,'FontSize',15)
% 
% %% Exit distribution of tagged particles if gathering source exit statistics
% if(~isempty(sink))
%     figure
%     frame_h = get(handle(gcf),'JavaFrame');
%     set(frame_h,'Maximized',1);
%     edges = [0.5:7.5]; % number of bins in histogram
%     h=histogram(output_matrix(:,5),edges);
%     xlabel('d (micron)')
%     ylabel('Number of particles')
%     E = h.BinEdges;
%     F = h.BinCounts;
%     xloc = E(1:end-1)+diff(E)/2;
%     text(xloc-0.5, F+5, string(F))
% end
% %% Plot of spectra
% % v1 = VideoWriter('particle_spectra.avi');
% % open(v1)
% % for i=1:numel(exit_mat(1,:))
% %     figure
% %     frame_h = get(handle(gcf),'JavaFrame');
% %     set(frame_h,'Maximized',1);
% % %     edges = [1 2 9 10 40 50 240 250 490 500]; % number of bins in histogram
% %     h=histogram(Np_inj(:,i).*(1-exit_mat(:,i))*1e7,edges);
% %     xlabel('d (micron) \times 10^{-1}')
% %     ylabel('Number of particles')
% % %     ylim([0 250])
% %     E = h.BinEdges;
% %     F = h.BinCounts;
% %     xloc = E(1:end-1)+diff(E)/2;
% %     text(xloc-0.5, F+5, string(F))
% %     frame=getframe(gcf); % leaving gcf out crops the frame in the movie.
% %     writeVideo(v1,frame)
% % end
% % close(v1)
% 
% end
% 
