
%% Proprietary code
%% Written by
%% Calina A. Copos (calinacopos@gmail.com | calinacopos@unc.edu)
%%
%% ``Stress fibers are embedded in a contractile cortical network" Nat. Material (2020)
%% Last updated: June 2020

close all;
clear;
clc;

%% Parameters
yloc = [-1 3];
xloc = [0 7];

Tmax = 2.5;
dt = 0.0001; 
drag = 0.025;

k1 = 1.0; gamma1 = 0.0;     % cytoskeleton parameters
k2 = 20; gamma2 = 2.0;      % force generator (stress fibers) parameters
k3 = k1; gamma3 = gamma1;   % fiber-cytoskeleton connectors parameters


%% Fresh configuration
[Top0,TopLinks,Bottom0,BottomLinks,Mid0,Tri] = formStructure;
Top = Top0;
Bottom = Bottom0;
Mid = Mid0;

%% Read previous configuration
AA = load('Configurations/relaxedconf_1_0.5_20_2_1_0.5.txt');

Nodes0 = [AA(:,1) AA(:,2)];
Nodes = [AA(:,3) AA(:,4)];
NT = 18; 
MT = length(Nodes0)-2*NT;
Top0 = Nodes0(1:NT,:);
Mid0 = Nodes0((NT+1):(NT+MT),:);
Bottom0 = Nodes0((NT+MT+1):end,:);   
Top = Nodes(1:NT,:);
Mid = Nodes((NT+1):(NT+MT),:);
Bottom = Nodes((NT+MT+1):end,:);
 
%% Compute initial time forces
Tcut = 0; Bcut = 0;

%FMid = elasticMidElements(Tri,Mid0,Mid,k1); % elastic only cytoskeleton
FMid = MidElements(Tri,Mid0,Mid,k1,gamma1);
FTop = contractileElements(Top0,Top,k2,gamma2,Tcut);
FBottom = contractileElements(Bottom0,Bottom,k2,gamma2,Bcut);

[FConnectTop,FConnectMidTop] = crossTopElements(Top0,Mid0,Top,Mid,TopLinks,k3,gamma3);
[FConnectBottom,FConnectMidBottom] = crossBottomElements(Bottom0,Mid0,Bottom,Mid,BottomLinks,k3,gamma3);

FMid = FMid + FConnectMidTop + FConnectMidBottom;
FTop = FTop + FConnectTop;
FBottom = FBottom + FConnectBottom; 

% Memory allocation
Ntotal = length(Top)+length(Mid)+length(Bottom);
NT = length(Top);
Nodes0 = [Top0;Mid0;Bottom0];
Nodes = [Top;Mid;Bottom];
FT = zeros(Ntotal,2);
FT = [FTop;FMid;FBottom];
FtopL = zeros(round(Tmax/dt),1);
FtopR = zeros(round(Tmax/dt),1);
FbottomL = zeros(round(Tmax/dt),1);
FbottomR = zeros(round(Tmax/dt),1);

%% Record initial traction stresses
% Only horizontal component of forces exert measured traction stresses
FtopL(1) = FT(1,1) + FT(NT+2,1) + FT(NT+8,1) + FT(NT+7,1);
FbottomL(1) = FT(NT+1,1) + FT(NT+5,1) + FT(NT+6,1) + FT(NT+length(Mid0)+1,1);
FtopR(1) = FT(NT,1) + FT(NT+4,1) + FT(NT+99,1) + FT(NT+97,1);
FbottomR(1) = FT(NT+3,1) + FT(NT+96,1) + FT(NT+95,1) + FT(Ntotal,1);

%% March through time
%Tcut = 1;
TopPos = zeros(round(Tmax/dt)-1,length(Top),2);
TopPos(1,:,:) = Top;
for i = 1:(Tmax/dt-1)
    Top = Nodes(1:NT,:);
    Mid = Nodes((NT+1):(NT+length(Mid)),:);
    Bottom = Nodes((NT+length(Mid)+1):end,:);
    % Contractile forces in top & bottom row 
    %if(i>5000) Tcut = 2; end
    %if(i>15000) Tcut = 2; end%Bcut = 3; end
    FTop = contractileElements(Top0,Top,k2,gamma2,Tcut);
    FBottom = contractileElements(Bottom0,Bottom,k2,gamma2,Bcut);
    %FMid = elasticMidElements(Tri,Mid0,Mid,k1); % elastic cytoskeleton
    FMid = MidElements(Tri,Mid0,Mid,k1,gamma1); % contractile cytoskeleton
    
    % Elastic forces between top and under top row
    % Perform ablation/shave only at the top
    %if(i>=1)
    %    [FConnectTop,FConnectMidTop] = crossTopElements(Top0,Mid0,Top,Mid,TopLinks,0.0,0.0);
    %else
        [FConnectTop,FConnectMidTop] = crossTopElements(Top0,Mid0,Top,Mid,TopLinks,k3,gamma3);
    %end
    [FConnectBottom,FConnectMidBottom] = crossBottomElements(Bottom0,Mid0,Bottom,Mid,BottomLinks,k3,gamma3);

    FMid = FMid + FConnectMidTop + FConnectMidBottom;
    FTop = FTop + FConnectTop;
    FBottom = FBottom + FConnectBottom; 
    
    % Total forces
    FT = zeros(Ntotal,2);
    FT = [FTop;FMid;FBottom];
    FTT(i,:,:) = FT;
    
    dXdt = zeros(Ntotal,2);
    dXdt = FT/drag; 
    
    % Adherent spots stay put
    dXdt(1,:) = [0 0];
    dXdt(NT+1,:) = [0 0];
    dXdt(NT+2,:) = [0 0];
    dXdt(NT+5,:) = [0 0];
    dXdt(NT+6,:) = [0 0];
    dXdt(NT+7,:) = [0 0];
    dXdt(NT+8,:) = [0 0];
    dXdt(NT+length(Mid0)+1,:) = [0 0];
    
    dXdt(NT,:) = [0 0];
    dXdt(NT+3,:) = [0 0];
    dXdt(NT+4,:) = [0 0];
    dXdt(NT+95,:) = [0 0];
    dXdt(NT+96,:) = [0 0];
    dXdt(NT+97,:) = [0 0];
    dXdt(NT+length(Mid0),:) = [0 0];
    dXdt(Ntotal,:) = [0 0];
    
    % update in time
    Nodes(:,1) = Nodes(:,1) + dXdt(:,1)*dt;
    Nodes(:,2) = Nodes(:,2) + dXdt(:,2)*dt;
    
    Top = Nodes(1:NT,:);
    Mid = Nodes((NT+1):(NT+length(Mid)),:);
    Bottom = Nodes((NT+length(Mid)+1):end,:);

    TopPos(i,:,:) = Top;
    % Plot
%    xlim([-0.5 7.5]); ylim([-1.5 3.5]); hold on;
%    for id=2:(NT-1)
%         iid=TopLinks(id,1);
%         line([Top(id,1),Mid(iid,1)],[Top(id,2),Mid(iid,2)],'color','b');
%         iid=TopLinks(id,2);
%         line([Top(id,1),Mid(iid,1)],[Top(id,2),Mid(iid,2)],'color','b');
%         
%         iiid=BottomLinks(id,1);
%         line([Bottom(id,1),Mid(iiid,1)],[Bottom(id,2),Mid(iiid,2)],'color','b');
%         iiid=BottomLinks(id,2);
%         line([Bottom(id,1),Mid(iiid,1)],[Bottom(id,2),Mid(iiid,2)],'color','b');
%     end
%     line([Top(1,1),Mid(2,1)],[Top(1,2),Mid(2,2)],'color','b');
%     line([Top(end,1),Mid(4,1)],[Top(end,2),Mid(4,2)],'color','b');
%     line([Bottom(1,1),Mid(1,1)],[Bottom(1,2),Mid(1,2)],'color','b');
%     line([Bottom(end,1),Mid(3,1)],[Bottom(end,2),Mid(3,2)],'color','b');
%     line(Top(:,1),Top(:,2),'color',[0 0 1],'linewidth',2);
%     line(Bottom(:,1),Bottom(:,2),'color',[0 0 1],'linewidth',2);
%     
%     scatter(Top(:,1),Top(:,2),100,'ob','fill','b'); 
%     scatter(Bottom(:,1),Bottom(:,2),100,'ob','fill','b');
%     scatter(Mid(:,1),Mid(:,2),100,'ok','fill','k');
%     triplot(Tri,Mid(:,1),Mid(:,2),'k');
%    set(gca,'fontsize',30); set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]); box on;
%    set(gca,'plotBoxAspectRatio',[8 5 1]); 
% 
% %     
%      quiver(Nodes(:,1),Nodes(:,2),FT(:,1),FT(:,2),'color',[0.75 0 0.75],'linewidth',5);
% %     %quiver(Top(:,1),Top(:,2),FConnectTop(:,1),FConnectTop(:,2),'linewidth',2);
% %     %quiver(UnderTop(:,1),UnderTop(:,2),FConnectUnderTop(:,1),FConnectUnderTop(:,2),'linewidth',2);
% %     %scatter(Nodes0(:,1),Nodes0(:,2),100,'ro','fill','r');
% %     %scatter(AA(:,3),AA(:,4),100,'ro','fill','r');
% %     %scatter(Nodes(:,1),Nodes(:,2),100,'kd','fill','k');
% %     set(gca,'fontsize',30); set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]); box on;
% %     set(gca,'plotBoxAspectRatio',[8 5 1]); 
% %     hold off; 
% % %     %pause; 
%       currFrame = getframe(gcf);
%      cla;

%    % Only horizontal component of forces exert measured traction stresses
%    % force computation with no middle strip
    FtopL(i+1) = FT(1,1) + FT(NT+2,1) + FT(NT+8,1) + FT(NT+7,1);
    FbottomL(i+1) = FT(NT+1,1) + FT(NT+5,1) + FT(NT+6,1) + FT(NT+length(Mid0)+1,1);
    FtopR(i+1) = FT(NT,1) + FT(NT+4,1) + FT(NT+99,1) + FT(NT+97,1);
    FbottomR(i+1) = FT(NT+3,1) + FT(NT+96,1) + FT(NT+95,1) + FT(Ntotal,1);
     
    if((abs(sum(FT(:,1)))>1e-13) || (abs(sum(FT(:,2)))>1e-13) ) 
        sprintf('TROUBLE!!!\f')
    end
end

%% Plot traction forces
Ftraction=FtopL+FbottomL+abs(FtopR)+abs(FbottomR);
abs(Ftraction(15000)-Ftraction(1))/Ftraction(1);
abs(Ftraction(end)-Ftraction(15000))/Ftraction(15000);
FTLoss = abs(Ftraction(end)-Ftraction(1));
% figure(2);
% plot(linspace(0,Tmax,length(Ftraction)),FtopL,'--','color','r','linewidth',2); hold on;
% plot(linspace(0,Tmax,length(Ftraction)),abs(FtopR),'-','color','k','linewidth',2);
% plot(linspace(0,Tmax,length(Ftraction)),FbottomL,'-','color','g','linewidth',2);hold on;
% plot(linspace(0,Tmax,length(Ftraction)),abs(FbottomR),'-.','color','b','linewidth',2);
% lgd = legend('top left','top right','bottom left','bottom right');
% set(gca,'fontsize',30);
% lgd.FontSize = 20; box on;
% xlabel('Time','interpreter','latex');
% ylabel('Force','interpreter','latex');
% ylim([0 5]);