function F = contractileElements(X0,X,k,gamma,cut)
    N = length(X);

    % setup neighbor list
    id = 1:N;
    left = id-1;
    left(1) = 1;
    right = id+1;
    right(end) = N;
    
    % setup edge list
    edgeList = zeros(N,2);
    for i=1:N
        edgeList(i,:) = [i-1 i];
    end
    edge_L = edgeList(:,1);
    edge_L(1) = 0;
    edge_R = edgeList(:,2);
    edge_R(end) = 0;
    
    GAMMA = gamma*ones(N,1);
    KK = k*ones(N,1);
    
    % compute reference and current distance to left and right neighbors
    ref_dist_L = sqrt((X0(:,1)-X0(left,1)).^2+(X0(:,2)-X0(left,2)).^2);
    ref_dist_R = sqrt((X0(right,1)-X0(:,1)).^2+(X0(right,2)-X0(:,2)).^2);
    def_dist_L = sqrt((X(:,1)-X(left,1)).^2+(X(:,2)-X(left,2)).^2);
    def_dist_R = sqrt((X(right,1)-X(:,1)).^2+(X(right,2)-X(:,2)).^2);
    
    % compute ds per structure node
    ds = 0.5*(def_dist_L+def_dist_R);
    
    % tangent vector for left and right neighbors
    tan_L(:,1) = -(X(:,1)-X(left,1))./def_dist_L;
    tan_L(:,2) = -(X(:,2)-X(left,2))./def_dist_L;
    tan_R(:,1) = (X(right,1)-X(:,1))./def_dist_R;
    tan_R(:,2) = (X(right,2)-X(:,2))./def_dist_R;
    
    % constitutive law: elastic, contractile material
    GAMMA_R = GAMMA-(edge_R==0); KK_R = KK-(edge_R==0);
    GAMMA_L = GAMMA-(edge_L==0); KK_L = KK-(edge_L==0);
    
    % make cut 
    if cut==1 % cut 1 performed between links 4-5
        GAMMA_R(4) = 0.02;%0.1;%0.2;%0.1*GAMMA_R(2);  % pill: 0.1
        GAMMA_L(5) = 0.02;%0.1;%0.2;%0.1*GAMMA_L(3);  % pill: 0.1
        
        KK_R(4) = 0.2;%0.2;%0.1*KK_R(2);        % pill: 0.0
        KK_L(5) = 0.2;%0.2;%0.1*KK_L(3);        % pill: 0.0
    end
    if cut==2 % cut 1 performed between links 4-5 followed by a second cut
        % on links 14-15
        GAMMA_R(4) = 0.02;%0.1*GAMMA_R(2);
        GAMMA_L(5) = 0.02;%0.1*GAMMA_L(3);
       
        KK_R(4) = 0.2;%0.1*KK_R(2);
        KK_L(5) = 0.2;%0.1*KK_L(3);
        
        GAMMA_R(14) = 0.02;%0.1*GAMMA_R(3);
        GAMMA_L(15) = 0.02;%0.1*GAMMA_L(4);
        
        KK_R(14) = 0.2;%0.1*KK_R(3);
        KK_L(15) = 0.2;%0.1*KK_L(4);
    end
    if cut==3
        GAMMA_R(14) = 0.02;%0.1*GAMMA_R(7);
        GAMMA_L(15) = 0.02;%0.1*GAMMA_L(8);
       
        KK_R(14) = 0.2;%0.1*KK_R(7);
        KK_L(15) = 0.2;%0.1*KK_L(8);
    end
    if cut==4
        pp = zeros(6,2);
        pp(1,:) = [7 8]; pp(2,:) = [5 6]; pp(3,:) = [3 4];
        pp(4,:) = [11 12]; pp(5,:) = [13 14]; pp(6,:) = [15 16];
        for kk=1:6
            GAMMA_R(pp(kk,1)) = 0.2;%0.1*GAMMA_R(2);
            GAMMA_L(pp(kk,2)) = 0.2;%0.1*GAMMA_L(3);

            KK_R(pp(kk,1)) = 0.2;%0.1*KK_R(2);
            KK_L(pp(kk,2)) = 0.2;%0.1*KK_L(3);
        end
    end
    if cut==5
        pp = zeros(7,2);
        pp(1,:) = [7 8]; pp(2,:) = [5 6]; pp(3,:) = [3 4];
        pp(4,:) = [11 12]; pp(5,:) = [13 14]; pp(6,:) = [15 16];
        pp(7,:) = [9 10];
        for kk=1:7
            GAMMA_R(pp(kk,1)) = 0.2;%0.1*GAMMA_R(2);
            GAMMA_L(pp(kk,2)) = 0.2;%0.1*GAMMA_L(3);

            KK_R(pp(kk,1)) = 0.2;%0.1*KK_R(2);
            KK_L(pp(kk,2)) = 0.2;%0.1*KK_L(3);
        end
    end
    if cut==6 %% mid cut performed between links 9-10
        GAMMA_R(9) = 0.02;%0.1*GAMMA_R(2);
        GAMMA_L(10) = 0.02;%0.1*GAMMA_L(3);
       
        KK_R(9) = 0.2;%0.1*KK_R(2);
        KK_L(10) = 0.2;%0.1*KK_L(3);
    end
    
   
    F(:,1) = (GAMMA_R+KK_R.*(def_dist_R./ref_dist_R - 1.0)).*tan_R(:,1);
    F(:,1) = F(:,1) + (GAMMA_L+KK_L.*(def_dist_L./ref_dist_L - 1.0)).*tan_L(:,1);
    %F(:,1) = F(:,1)./ds;

    F(:,2) = (GAMMA_R+KK_R.*(def_dist_R./ref_dist_R - 1.0)).*tan_R(:,2);
    F(:,2) = F(:,2) + (GAMMA_L+KK_L.*(def_dist_L./ref_dist_L - 1.0)).*tan_L(:,2);
    %F(:,2) = F(:,2)./ds;
    
    F(1,:) = (GAMMA_R(1)+KK_R(1).*(def_dist_R(1)./ref_dist_R(1) - 1.0)).*tan_R(1,:);
    F(end,:) = (GAMMA_L(end)+KK_L(end).*(def_dist_L(end)./ref_dist_L(end) - 1.0)).*tan_L(end,:);
end

