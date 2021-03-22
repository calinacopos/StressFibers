function F = MidElements(T,X0,X,kk,gamma)
    
    N = length(X);
    Ntris = length(T);
    
    % zero vectors
    Fc = zeros(N,2);
    dA_T = zeros(Ntris,1);
    dA = zeros(N,1);
    
    % fetch edges
    TR = TriRep(T,X0(:,1),X0(:,2));
    E = edges(TR);
    
    % area elements
    dA_T = 0.5*( (X(T(:,2),1)-X(T(:,1),1)).*(X(T(:,3),2)-X(T(:,1),2)) - ...
        (X(T(:,2),2)-X(T(:,1),2)).*(X(T(:,3),1)-X(T(:,1),1)) );
    for k=1:Ntris
        dA(T(k,1)) = dA(T(k,1)) + dA_T(k)/3;
        dA(T(k,2)) = dA(T(k,2)) + dA_T(k)/3;
        dA(T(k,3)) = dA(T(k,3)) + dA_T(k)/3; 
    end
    
    % for each edge
    for k=1:length(E)
        id1 = E(k,1);
        id2 = E(k,2);
        dl0 = sqrt((X0(id2,1)-X0(id1,1))^2 + (X0(id2,2)-X0(id1,2))^2);
        dl = sqrt((X(id2,1)-X(id1,1))^2 + (X(id2,2)-X(id1,2))^2);
        dk = (8*kk/(3*dl0))*(dA(id1)+dA(id2))/2;
        
        ddA = (dA(id1)+dA(id2))/2;
    
        Fc(id1,1) = Fc(id1,1) + ( gamma + (dk/ddA)*(dl/dl0 - 1.0) )*(X(id2,1)-X(id1,1))/dl;
        Fc(id1,2) = Fc(id1,2) + ( gamma + (dk/ddA)*(dl/dl0 - 1.0) )*(X(id2,2)-X(id1,2))/dl;
    
        Fc(id2,1) = Fc(id2,1) + ( gamma + (dk/ddA)*(dl/dl0 - 1.0) )*(X(id1,1)-X(id2,1))/dl;
        Fc(id2,2) = Fc(id2,2) + ( gamma + (dk/ddA)*(dl/dl0 - 1.0) )*(X(id1,2)-X(id2,2))/dl;
    end
        
    F(:,1) = Fc(:,1); 
    F(:,2) = Fc(:,2);
end
