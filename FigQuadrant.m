% Figure for quadrant distribution
figure(5); hold on;
FTLoss = abs(Ftraction(end)-Ftraction(1));
TL = abs(FtopL(end)-FtopL(1))/FTLoss;
TR = abs(FtopR(end)-FtopR(1))/FTLoss;
BL = abs(FbottomL(end)-FbottomL(1))/FTLoss;
BR = abs(FbottomR(end)-FbottomR(1))/FTLoss;

subplot(1,2,1); hold on;
scatter(1,TL*100,5000,'o','MarkerEdgeColor',[191 191 0]/255,'MarkerFaceColor',[191 191 0]/255);
plot(linspace(0.5,1.5,100),TL*100*ones(100,1),'linewidth',2,'color',[80 80 80]/255)
scatter(2,TR*100,5000,'s','MarkerEdgeColor',[0 0 255]/255,'MarkerFaceColor',[0 0 255]/255);
plot(linspace(1.5,2.5,100),TR*100*ones(100,1),'linewidth',2,'color',[80 80 80]/255)
ylim([0 80]); xlim([0 2.5]);
set(gca,'XTickLabel',[]);
set(gca,'linewidth',4); box off;
set(gca,'fontsize',30);

subplot(1,2,2); hold on;
scatter(1,(TL+TR)*100,5000,'o','MarkerEdgeColor',[0 128 0]/255,'MarkerFaceColor',[0 128 0]/255);
plot(linspace(0.5,1.5,100),(TL+TR)*100*ones(100,1),'linewidth',2,'color',[80 80 80]/255)
scatter(2,(BL+BR)*100,5000,'d','MarkerEdgeColor',[191 0 191]/255,'MarkerFaceColor',[191 0 191]/255);
plot(linspace(1.5,2.5,100),(BL+BR)*100*ones(100,1),'linewidth',2,'color',[80 80 80]/255)
ylim([0 100]); xlim([0 2.5]);
set(gca,'XTickLabel',[]);
set(gca,'linewidth',4); box off;
set(gca,'fontsize',30);
 