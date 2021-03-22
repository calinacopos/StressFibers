function [Top0,TopLinks,Bottom0,BottomLinks,Mid0,Tri] = formStructure
    p = load('Mesh/midconf_PointsN99.txt');
    t = load('Mesh/midconf_TrianglesN99.txt');
    bdd = load('Mesh/midconf_BoundaryN99.txt');
    
    x = p(:,1); y = p(:,2);
    Mid0 = [x y];
    Tri = t;
    
    topId = find(y>1.9);
    topX = x(topId);
    topX = sort(topX);
    N = length(topId)+1;
    Top0 = zeros(N,2);
    Top0(:,2) = 3;
    Top0(1,1) = 0; Top0(end,1) = 7;    
    for i=1:(length(topId)-1)
        Top0(i+1,1) = (topX(i)+topX(i+1))/2;
    end
    
    % TopLinks (organized configuration)
    TopLinks = zeros(N,2);
    TopLinks(1,:) = [2 -1];
    TopLinks(2,:) = [2 9];
    TopLinks(3,:) = [9 15];
    TopLinks(4,:) = [15 21];
    TopLinks(5,:) = [21 27];
    TopLinks(6,:) = [27 33];
    TopLinks(7,:) = [33 39];
    TopLinks(8,:) = [39 45];
    TopLinks(9,:) = [45 51];
    TopLinks(10,:) = [51 57];
    TopLinks(11,:) = [57 63];
    TopLinks(12,:) = [63 69];
    TopLinks(13,:) = [69 75];
    TopLinks(14,:) = [75 81];
    TopLinks(15,:) = [81 87];
    TopLinks(16,:) = [87 93];
    TopLinks(17,:) = [93 4];
    TopLinks(end,:) = [4 -1];
    
    botId = find(y<0.1);
    botX = x(botId);
    botX = sort(botX);
    N = length(botId)+1;
    Bottom0 = zeros(N,2);
    Bottom0(:,2) = -1;
    Bottom0(1,1) = 0; Bottom0(end,1) = 7;    
    for i=1:(length(botId)-1)
        Bottom0(i+1,1) = (botX(i)+botX(i+1))/2;
    end

    % BottomLinks (organized configuration)
    BottomLinks = zeros(N,2);
    BottomLinks(1,:) = [1 -1];
    BottomLinks(2,:) = [1 10];
    BottomLinks(3,:) = [10 16];
    BottomLinks(4,:) = [16 22];
    BottomLinks(5,:) = [22 28];
    BottomLinks(6,:) = [28 34];
    BottomLinks(7,:) = [34 40];
    BottomLinks(8,:) = [40 46];
    BottomLinks(9,:) = [46 52];
    BottomLinks(10,:) = [52 58];
    BottomLinks(11,:) = [58 64];
    BottomLinks(12,:) = [64 70];
    BottomLinks(13,:) = [70 76];
    BottomLinks(14,:) = [76 82];
    BottomLinks(15,:) = [82 88];
    BottomLinks(16,:) = [88 94];
    BottomLinks(17,:) = [94 3];
    BottomLinks(end,:) = [3 -1];
end
