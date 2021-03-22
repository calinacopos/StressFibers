% Figure for configuration
figure(2); 

subplot(1,2,2);
plot(linspace(0,Tmax,length(Ftraction)),abs(Ftraction),'-k','linewidth',10);
xlim([0,Tmax]); ylim([10 16]);
set(gca,'fontsize',30); set(gca,'linewidth',4); box off; hold on;

% cut1 = abs(Ftraction(10000)-Ftraction(1))/Ftraction(1);
% cut2 = abs(Ftraction(end)-Ftraction(10000))/Ftraction(1);
% set(gca,'fontsize',30); set(gca,'linewidth',4); box off; hold on;
% scatter(1,cut1*100,5000,'o','MarkerEdgeColor',[255 0 0]/255,'MarkerFaceColor',[255 0 0]/255);
% plot(linspace(0.5,1.5,100),cut1*100*ones(100,1),'linewidth',2,'color',[80 80 80]/255)
% scatter(2,cut2*100,5000,'o','MarkerEdgeColor',[0 0 255]/255,'MarkerFaceColor',[0 0 255]/255);
% plot(linspace(1.5,2.5,100),cut2*100*ones(100,1),'linewidth',2,'color',[80 80 80]/255);
% xlim([0 2.5]); ylim([0 50]); set(gca,'XTickLabel',[]);

subplot(1,2,1);
xlim([-0.5 7.5]); ylim([-1.5 3.5]); hold on;
for id=2:(NT-1)
    iid=TopLinks(id,1);
    line([Top(id,1),Mid(iid,1)],[Top(id,2),Mid(iid,2)],'color','b');
    iid=TopLinks(id,2);
    line([Top(id,1),Mid(iid,1)],[Top(id,2),Mid(iid,2)],'color','b');

    iiid=BottomLinks(id,1);
    line([Bottom(id,1),Mid(iiid,1)],[Bottom(id,2),Mid(iiid,2)],'color','b');
    iiid=BottomLinks(id,2);
    line([Bottom(id,1),Mid(iiid,1)],[Bottom(id,2),Mid(iiid,2)],'color','b');
end
line([Top(1,1),Mid(2,1)],[Top(1,2),Mid(2,2)],'color','b');
line([Top(end,1),Mid(4,1)],[Top(end,2),Mid(4,2)],'color','b');
line([Bottom(1,1),Mid(1,1)],[Bottom(1,2),Mid(1,2)],'color','b');
line([Bottom(end,1),Mid(3,1)],[Bottom(end,2),Mid(3,2)],'color','b');
line(Top(:,1),Top(:,2),'color',[0 0 1],'linewidth',2);
line(Bottom(:,1),Bottom(:,2),'color',[0 0 1],'linewidth',2);

scatter(Top(:,1),Top(:,2),100,'o','markerfacecolor',[255 255 255]/255,'markeredgecolor','b'); 
scatter(Bottom(:,1),Bottom(:,2),100,'o','markerfacecolor',[255 255 255]/255,'markeredgecolor','b');
scatter(Mid(:,1),Mid(:,2),100,'o','markerfacecolor',[255 255 255]/255,'markeredgecolor','k');
triplot(Tri,Mid(:,1),Mid(:,2),'k');
set(gca,'fontsize',30); set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]); box on;
set(gca,'plotBoxAspectRatio',[8 5 1]); 

%mm = find(Mid(:,2)>0.7 & Mid(:,2)<1.3);
%scatter(Mid(mm,1),Mid(mm,2),100,'ok','fill','k');
scatter(Top(1,1),Top(1,2),100,'ob','fill','b');
scatter(Top(NT,1),Top(NT,2),100,'ob','fill','b');
scatter(Bottom(NT,1),Bottom(NT,2),100,'ob','fill','b');
scatter(Bottom(1,1),Bottom(1,2),100,'ob','fill','b');
scatter(Nodes(NT+3,1),Nodes(NT+3,2),100,'ko','fill','k');
scatter(Nodes(NT+95,1),Nodes(NT+95,2),100,'ko','fill','k');
scatter(Nodes(NT+96,1),Nodes(NT+96,2),100,'ko','fill','k');
scatter(Nodes(NT+97,1),Nodes(NT+97,2),100,'ko','fill','k');
scatter(Nodes(NT+99,1),Nodes(NT+99,2),100,'ko','fill','k');
scatter(Nodes(NT+4,1),Nodes(NT+4,2),100,'ko','fill','k');
scatter(Nodes(NT+1,1),Nodes(NT+1,2),100,'ko','fill','k');
scatter(Nodes(NT+7,1),Nodes(NT+7,2),100,'ko','fill','k');
scatter(Nodes(NT+8,1),Nodes(NT+8,2),100,'ko','fill','k');
scatter(Nodes(NT+2,1),Nodes(NT+2,2),100,'ko','fill','k');
scatter(Nodes(NT+5,1),Nodes(NT+5,2),100,'ko','fill','k');
scatter(Nodes(NT+6,1),Nodes(NT+6,2),100,'ko','fill','k');