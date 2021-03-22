% Figure for distance of retraction away from the centered cut
figure(5); hold on;

Top_rescaled_to_experiments = Top/7*37; % experiments report 37um for pill 
                                        % shape, while I use 7 a.u.
Top0_rescaled_to_experiments = Top0/7*37;
                                        
for i=1:18
    distRetract(i)  = abs( TopPos(end,i,1)-TopPos(1,i,1) );
    distCut(i)      = abs( Top_rescaled_to_experiments(i,1)-0.5*(Top0_rescaled_to_experiments(9,1)+Top0_rescaled_to_experiments(10,1)) );
end
scatter(distCut,distRetract,1000,'o','MarkerEdgeColor','k','MarkerFaceColor',[204 204 204]/255,'linewidth',2);
set(gca,'fontsize',30); set(gca,'linewidth',4); box off; hold on;