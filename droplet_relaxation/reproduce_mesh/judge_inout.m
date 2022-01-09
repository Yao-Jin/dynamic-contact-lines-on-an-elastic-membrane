function [inout] = judge_inout(node,elem,bd_droplet)

    % judge if an element inside or outside the droplet bounded by bd_droplet
    
    node_inout = inpolygon(node(:,1),node(:,2),bd_droplet(:,1),bd_droplet(:,2));
    inout = node_inout(elem(:,1)).*node_inout(elem(:,2)).*node_inout(elem(:,3));
    
end
   