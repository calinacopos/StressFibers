function [FTop,FMidTop] = crossTopElements(Top0,Mid0,Top,Mid,TopLinks,k,gamma);
    % zero forces
    N = length(Top0);
    FTop = zeros(length(Top0),2);
    FMidTop = zeros(length(Mid0),2);
    
    % for interior points
    id = 2:N-1;
    
    % loop over pair in links
    for j=1:2
        F = zeros(length(id),2);
        v = zeros(length(id),2);
        
        iid = TopLinks(id,j);
    
        ref_dist = sqrt((Top0(id,1)-Mid0(iid,1)).^2+(Top0(id,2)-Mid0(iid,2)).^2);
        def_dist = sqrt((Top(id,1)-Mid(iid,1)).^2+(Top(id,2)-Mid(iid,2)).^2);
    
        % directionality vector on top row
        v(id,1) = (Mid(iid,1)-Top(id,1))./def_dist;
        v(id,2) = (Mid(iid,2)-Top(id,2))./def_dist;
    
        % force vector field on top row
        F(:,1) = (gamma+k.*(def_dist./ref_dist - 1.0)).*v(id,1);
        F(:,2) = (gamma+k.*(def_dist./ref_dist - 1.0)).*v(id,2);
        FTop(id,1) = FTop(id,1) + F(:,1);
        FTop(id,2) = FTop(id,2) + F(:,2);
    
        FMidTop(iid,1) = FMidTop(iid,1) - F(:,1);
        FMidTop(iid,2) = FMidTop(iid,2) - F(:,2);
    end
    
    % special cases: endpoints
    id = 1; iid = 2;
    ref_dist = sqrt((Top0(id,1)-Mid0(iid,1)).^2+(Top0(id,2)-Mid0(iid,2)).^2);
    def_dist = sqrt((Top(id,1)-Mid(iid,1)).^2+(Top(id,2)-Mid(iid,2)).^2);

    % directionality vector on top row
    v(id,1) = (Mid(iid,1)-Top(id,1))./def_dist;
    v(id,2) = (Mid(iid,2)-Top(id,2))./def_dist;

    % force vector field on top row
    F = zeros(2,1);
    F(1) = (gamma+k.*(def_dist./ref_dist - 1.0)).*v(id,1);
    F(2) = (gamma+k.*(def_dist./ref_dist - 1.0)).*v(id,2);
    FTop(id,1) = FTop(id,1) + F(1);
    FTop(id,2) = FTop(id,2) + F(2);
    
    FMidTop(iid,1) = FMidTop(iid,1) - F(1);
    FMidTop(iid,2) = FMidTop(iid,2) - F(2);
    
    id = N; iid = 4;
    ref_dist = sqrt((Top0(id,1)-Mid0(iid,1)).^2+(Top0(id,2)-Mid0(iid,2)).^2);
    def_dist = sqrt((Top(id,1)-Mid(iid,1)).^2+(Top(id,2)-Mid(iid,2)).^2);

    % directionality vector on top row
    v(id,1) = (Mid(iid,1)-Top(id,1))./def_dist;
    v(id,2) = (Mid(iid,2)-Top(id,2))./def_dist;

    % force vector field on top row
    F = zeros(2,1);
    F(1) = (gamma+k.*(def_dist./ref_dist - 1.0)).*v(id,1);
    F(2) = (gamma+k.*(def_dist./ref_dist - 1.0)).*v(id,2);
    FTop(id,1) = FTop(id,1) + F(1);
    FTop(id,2) = FTop(id,2) + F(2);

    FMidTop(iid,1) = FMidTop(iid,1) - F(1);
    FMidTop(iid,2) = FMidTop(iid,2) - F(2);
end

