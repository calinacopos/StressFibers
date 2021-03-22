% Figure for retraction away from the centered cut

set(gca,'Color','k'); hold on;
for i=1:18
plot(linspace(0,Tmax,length(TopPos)),TopPos(:,i,1),'-w','linewidth',4);
end
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
xlim([0.4 1]);