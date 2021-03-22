function [FBottom,FMidBottom] = crossBottomElements(Bottom0,Mid0,Bottom,Mid,BottomLinks,k,gamma);
    % zero forces
    N = length(Bottom0);
    FBottom = zeros(length(Bottom0),2);
    FMidBottom = zeros(length(Mid0),2);
    
    % for interior points
    id = 2:N-1;
    
    % loop over pair in links
    for j=1:2
        F = zeros(length(id),2);
        v = zeros(length(id),2);
        
        iid = BottomLinks(id,j);
    
        ref_dist = sqrt((Bottom0(id,1)-Mid0(iid,1)).^2+(Bottom0(id,2)-Mid0(iid,2)).^2);
        def_dist = sqrt((Bottom(id,1)-Mid(iid,1)).^2+(Bottom(id,2)-Mid(iid,2)).^2);
    
        % directionality vector on top row
        v(id,1) = (Mid(iid,1)-Bottom(id,1))./def_dist;
        v(id,2) = (Mid(iid,2)-Bottom(id,2))./def_dist;
    
        % force vector field on top row
        F(:,1) = (gamma+k.*(def_dist./ref_dist - 1.0)).*v(id,1);
        F(:,2) = (gamma+k.*(def_dist./ref_dist - 1.0)).*v(id,2);
        FBottom(id,1) = FBottom(id,1) + F(:,1);
        FBottom(id,2) = FBottom(id,2) + F(:,2);
    
        FMidBottom(iid,1) = FMidBottom(iid,1) - F(:,1);
        FMidBottom(iid,2) = FMidBottom(iid,2) - F(:,2);
    end
    
    % special cases: endpoints
    id = 1; iid = 1;
    ref_dist = sqrt((Bottom0(id,1)-Mid0(iid,1)).^2+(Bottom0(id,2)-Mid0(iid,2)).^2);
    def_dist = sqrt((Bottom(id,1)-Mid(iid,1)).^2+(Bottom(id,2)-Mid(iid,2)).^2);

    % directionality vector on top row
    v(id,1) = (Mid(iid,1)-Bottom(id,1))./def_dist;
    v(id,2) = (Mid(iid,2)-Bottom(id,2))./def_dist;

    % force vector field on top row
    F = zeros(2,1);
    F(1) = (gamma+k.*(def_dist./ref_dist - 1.0)).*v(id,1);
    F(2) = (gamma+k.*(def_dist./ref_dist - 1.0)).*v(id,2);
    FBottom(id,1) = FBottom(id,1) + F(1);
    FBottom(id,2) = FBottom(id,2) + F(2);

    FMidBottom(iid,1) = FMidBottom(iid,1) - F(1);
    FMidBottom(iid,2) = FMidBottom(iid,2) - F(2);
    
    id = N; iid = 3;
    ref_dist = sqrt((Bottom0(id,1)-Mid0(iid,1)).^2+(Bottom0(id,2)-Mid0(iid,2)).^2);
    def_dist = sqrt((Bottom(id,1)-Mid(iid,1)).^2+(Bottom(id,2)-Mid(iid,2)).^2);

    % directionality vector on top row
    v(id,1) = (Mid(iid,1)-Bottom(id,1))./def_dist;
    v(id,2) = (Mid(iid,2)-Bottom(id,2))./def_dist;

    % force vector field on top row
    F = zeros(2,1);
    F(1) = (gamma+k.*(def_dist./ref_dist - 1.0)).*v(id,1);
    F(2) = (gamma+k.*(def_dist./ref_dist - 1.0)).*v(id,2);
    FBottom(id,1) = FBottom(id,1) + F(1);
    FBottom(id,2) = FBottom(id,2) + F(2);

    FMidBottom(iid,1) = FMidBottom(iid,1) - F(1);
    FMidBottom(iid,2) = FMidBottom(iid,2) - F(2);
end

