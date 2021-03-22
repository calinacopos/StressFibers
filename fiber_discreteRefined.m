close all;
clear;
clc;
%% Parameters
yloc = [-1 3];
xloc = [0 7];
% structure params for dumbbell
k1 = 1.0; gamma1 = 0.5;     %k1 = 0.044; gamma1 = 0.2;   % cytoskeleton params
k2 = 15; gamma2 = 2.0;      % force generator (stress fibers) params
k3 = k1; gamma3 = gamma1;   % fiber-cyto connectors params
%k3 = 2.0; gamma3 = 0.75;     % fiber-cyto connectors params
% % % structure params for pill
% k1 = 0.5; gamma1 = 0.2;
% k2 = 5.0; gamma2 = 1.0;
% k3 = 0.5; gamma3 = 0.2;
% temporal 
Tmax = 2.5;
dt = 0.0001; 
drag = 0.025;

%% Fresh configuration
[Top0,TopLinks,Bottom0,BottomLinks,Mid0,Tri] = formStructure;
Top = Top0;
Bottom = Bottom0;
Mid = Mid0;

%% Read previous configuration
%AA = load('Configurations/relaxedconf_0.044_20_2_1_0.5.txt');
%AA = load('Configurations/relaxedconf_0.5_0.2_20_2_1_0.5.txt');
%AA = load('Configurations/relaxedconf_midstrip_0.5_0.0_20_2_1_0.5.txt');
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
% N.B. Aug 20 (need to include interior points too)
FtopL(1) = FT(1,1) + FT(NT+2,1) + FT(NT+8,1) + FT(NT+7,1);
FbottomL(1) = FT(NT+1,1) + FT(NT+5,1) + FT(NT+6,1) + FT(NT+length(Mid0)+1,1);
FtopR(1) = FT(NT,1) + FT(NT+4,1) + FT(NT+99,1) + FT(NT+97,1);
FbottomR(1) = FT(NT+3,1) + FT(NT+96,1) + FT(NT+95,1) + FT(Ntotal,1);

% force computation with middle strip (5 points)
% FtopL(1) = FT(1,1) + FT(NT+2,1) + FT(NT+8,1) + FT(NT+7,1) + FT(NT+25,1) + FT(NT+43,1);
% FbottomL(1) = FT(NT+1,1) + FT(NT+5,1) + FT(NT+6,1) + FT(NT+length(Mid0)+1,1) + FT(NT+12,1);
% FtopR(1) = FT(NT,1) + FT(NT+4,1) + FT(NT+99,1) + FT(NT+97,1) + FT(NT+61,1);
% FbottomR(1) = FT(NT+3,1) + FT(NT+96,1) + FT(NT+95,1) + FT(Ntotal,1) + FT(NT+84,1);

% force computation with middle strip (10 points)
% FtopL(1) = FT(1,1) + FT(NT+2,1) + FT(NT+8,1) + FT(NT+7,1) + FT(NT+25,1) + FT(NT+37,1) + FT(NT+43,1) ;
% FbottomL(1) = FT(NT+1,1) + FT(NT+5,1) + FT(NT+6,1) + FT(NT+length(Mid0)+1,1) + FT(NT+12,1) + FT(NT+48,1);
% FtopR(1) = FT(NT,1) + FT(NT+4,1) + FT(NT+99,1) + FT(NT+97,1) + FT(NT+73,1) + FT(NT+79,1) + FT(NT+61,1);
% FbottomR(1) = FT(NT+3,1) + FT(NT+96,1) + FT(NT+95,1) + FT(Ntotal,1) + FT(NT+54,1) + FT(NT+84,1);

% Pill setup
%kk = randperm(length(Mid));
%kk = kk(1:50);
%id = find(kk==2|kk==3|kk==4|kk==5|kk==6);
%kk(id) = 1;

% Adhesive strip through `middle'
%mm = find(Mid(:,2)>0.5 & Mid(:,2)<1.4);
%rmm = randperm(length(mm),10);
%rmid = mm(rmm);
%rmid = mm;
%rmid = [54    61    43    48    84    73     12    25    79    37]; % 10 ponts
%rmid = [61    12    25    84    43];                                % 5 points

% Extra fibers
% xxx = zeros(10,2);
% xxx(1,:) = [93 82]; xxx(2,:) = [63 71]; xxx(3,:) = [56 40]; xxx(4,:) = [32 29]; 
% xxx(5,:) = [21 10]; xxx(6,:) = [81 88]; xxx(7,:) = [65 74]; xxx(8,:) = [57 58]; 
% xxx(9,:) = [45 35]; xxx(10,:) = [9 28]; 

% Extra organized fibers
%xxx = formOrganizedXXXFibers;

%% March through time
%figure(1);
%Tcut = 1;
TopPos = zeros(round(Tmax/dt)-1,length(Top),2);
TopPos(1,:,:) = Top;
for i = 1:(Tmax/dt-1)
    Top = Nodes(1:NT,:);
    Mid = Nodes((NT+1):(NT+length(Mid)),:);
    Bottom = Nodes((NT+length(Mid)+1):end,:);
    % Contractile forces in top & bottom row 
    if(i>5000) Tcut = 1; end
    %if(i>15000) Tcut = 2; end%Bcut = 3; end
    FTop = contractileElements(Top0,Top,k2,gamma2,Tcut);
    FBottom = contractileElements(Bottom0,Bottom,k2,gamma2,Bcut);
    %FMid = elasticMidElements(Tri,Mid0,Mid,k1); % elastic cytoskeleton
    FMid = MidElements(Tri,Mid0,Mid,k1,gamma1); % contractile cytoskeleton
    %FExtraFibers = extraCytoElements(xxx,Mid0,Mid,k1,gamma1); 
    
    % Elastic forces between top and under top row
    % Perform ablation/shave only at the top
%     if(i>=1)
%         [FConnectTop,FConnectMidTop] = crossTopElements(Top0,Mid0,Top,Mid,TopLinks,0.0,0.0);
%     else
        [FConnectTop,FConnectMidTop] = crossTopElements(Top0,Mid0,Top,Mid,TopLinks,k3,gamma3);
%    end
    [FConnectBottom,FConnectMidBottom] = crossBottomElements(Bottom0,Mid0,Bottom,Mid,BottomLinks,k3,gamma3);

    %FMid = FMid + FConnectMidTop + FConnectMidBottom + FExtraFibers;
    FMid = FMid + FConnectMidTop + FConnectMidBottom;
    FTop = FTop + FConnectTop;
    FBottom = FBottom + FConnectBottom; 
    
    % Total forces
    FT = zeros(Ntotal,2);
    FT = [FTop;FMid;FBottom];
    FTT(i,:,:) = FT;
    
    dXdt = zeros(Ntotal,2);
    dXdt = FT/drag; 
    
%    % Pill setup
%     %dXdt(NT+kk',:) = zeros(length(kk),2);
%     dXdt([1 7:NT],:) = zeros(length([1 7:NT]),2);
%     dXdt((NT+length(Mid)+1):end,:) = zeros(NT,2);
    
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
    
    % Adhesive strip through middle of the cytoplasmic network
    %dXdt(NT+rmid',:) = zeros(length(rmid),2);
         
    % thin cytoplasmic network
%      dXdt(NT+9,:) = [0 0];
%      dXdt(NT+21,:) = [0 0];
%      dXdt(NT+39,:) = [0 0];
%      dXdt(NT+51,:) = [0 0];
%      dXdt(NT+75,:) = [0 0];
%      dXdt(NT+81,:) = [0 0];

    % large cytoplasmic network
%     dXdt(NT+10,:) = [0 0];
%     dXdt(NT+16,:) = [0 0];
%     dXdt(NT+40,:) = [0 0];
%     dXdt(NT+70,:) = [0 0];
%     dXdt(NT+76,:) = [0 0];
%     dXdt(NT+94,:) = [0 0];
    
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
    
    % force computation with middle strip (5 points)
%     FtopL(i+1) = FT(1,1) + FT(NT+2,1) + FT(NT+8,1) + FT(NT+7,1) + FT(NT+25,1) + FT(NT+43,1);
%     FbottomL(i+1) = FT(NT+1,1) + FT(NT+5,1) + FT(NT+6,1) + FT(NT+length(Mid0)+1,1) + FT(NT+12,1);
%     FtopR(i+1) = FT(NT,1) + FT(NT+4,1) + FT(NT+99,1) + FT(NT+97,1) + FT(NT+61,1);
%     FbottomR(i+1) = FT(NT+3,1) + FT(NT+96,1) + FT(NT+95,1) + FT(Ntotal,1) + FT(NT+84,1);

    % Only horizontal component of forces exert measured traction stresses
    % force computation with middle strip (10 points)
%     FtopL(i+1) = FT(1,1) + FT(NT+2,1) + FT(NT+8,1) + FT(NT+7,1) + FT(NT+25,1) + FT(NT+37,1) + FT(NT+43,1);
%     FbottomL(i+1) = FT(NT+1,1) + FT(NT+5,1) + FT(NT+6,1) + FT(NT+length(Mid0)+1,1) + FT(NT+12,1) + FT(NT+48,1);
%     FtopR(i+1) = FT(NT,1) + FT(NT+4,1) + FT(NT+99,1) + FT(NT+97,1) + FT(NT+73,1) + FT(NT+79,1) +  FT(NT+61,1);
%     FbottomR(i+1) = FT(NT+3,1) + FT(NT+96,1) + FT(NT+95,1) + FT(Ntotal,1) + FT(NT+54,1) + FT(NT+84,1);
    
    
    if((abs(sum(FT(:,1)))>1e-13) || (abs(sum(FT(:,2)))>1e-13) ) 
        sprintf('TROUBLE!!!\f')
    end
    %sprintf('sum forces x-dir %E sum forces y-dir %E\f',sum(FT(:,1)),sum(FT(:,2)))
end

%% Plot traction forces
Ftraction=FtopL+FbottomL+abs(FtopR)+abs(FbottomR);
abs(Ftraction(end)-Ftraction(1))/Ftraction(1)
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