function [geom] = adaptLRbd(geom,para,sln)

% before reproduce the numerical mesh, redistribute the markers on the left and right boundaries
% redistribution is to guarantee the same ratio of distances between markers to the total length

    Nmem = para.Nmem; 
    lbd2node = geom.lbd2node; rbd2node = geom.rbd2node;
    leftbdy =[sln.pmemold(1,2); geom.node(lbd2node,2);1]; 
    rightbdy = [sln.pmemold(Nmem,2); geom.node(rbd2node,2);1]; 
    ratio = (leftbdy(2:end)-leftbdy(1:end-1))/(1-leftbdy(1));
    tmp = 1-sln.pmem(1,2);
    n = size(leftbdy,1);
    newleftbdy = zeros(n,1);
    newleftbdy(1) = sln.pmem(1,2);
    for i = 2:n
        newleftbdy(i) = newleftbdy(i-1)+ratio(i-1)*tmp;
    end
    geom.node(lbd2node,2) = newleftbdy(2:end-1);
    
    ratio = (rightbdy(2:end)-rightbdy(1:end-1))/(1-rightbdy(1));
    tmp = 1-sln.pmem(Nmem,2);
    n = size(rightbdy,1);
    newrightbdy = zeros(n,1);
    newrightbdy(1) = sln.pmem(Nmem,2);
    for i = 2:n
        newrightbdy(i) = newrightbdy(i-1)+ratio(i-1)*tmp;
    end
    geom.node(rbd2node,2) = newrightbdy(2:end-1);    
    
    plot(ones(size(rightbdy)), rightbdy, '-or');
    hold on;
    plot(ones(size(newrightbdy)), newrightbdy, '-sb');
    hold off;
    
end