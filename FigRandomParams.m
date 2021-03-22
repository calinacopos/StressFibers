scatter(rand(N,1),MultipleCutLoss(:,4),1000,'MarkerEdgeColor','r','MarkerFaceColor',[0.7 0.7 0.7],'linewidth',1)
box on;
hold on;
scatter(rand(N,1)+3,MultipleCutLoss(:,5),1000,'MarkerEdgeColor','b','MarkerFaceColor',[0.7 0.7 0.7],'linewidth',1)
set(gcf,'color','w');
set(gca,'fontsize',20,'fontname','times');
set(gca,'xtick',[]);
ylim([0 1]); xlim([-5 8]);
plot(linspace(-0.5,1.5,10),mean(MultipleCutLoss(:,4))*ones(10,1),'-r','linewidth',5);
plot(linspace(2.5,4.5,10),mean(MultipleCutLoss(:,5))*ones(10,1),'-b','linewidth',5);


save ../Results_Rebuttal/Data/Fig4f_RandomElasticities.mat MultipleCutLoss