function kappa = getcurvature(p)

    N = size(p,1);
    len = sqrt((p(2:end,1)-p(1:end-1,1)).^2+(p(2:end,2)-p(1:end-1,2)).^2);
    
    tauelem = zeros(N-1,2);
    tauelem(:,1) = (p(2:end,1)-p(1:end-1,1))./len;
    tauelem(:,2) = (p(2:end,2)-p(1:end-1,2))./len;
    
    taunode = zeros(N,2);
    taunode(1,1:2) = tauelem(1,1:2);
    taunode(2:N-1,1) = len(1:N-2)./(len(1:N-2)+len(2:N-1)).*tauelem(1:N-2,1)...
                       + len(2:N-1)./(len(1:N-2)+len(2:N-1)).*tauelem(2:N-1,1);
    taunode(2:N-1,2) = len(1:N-2)./(len(1:N-2)+len(2:N-1)).*tauelem(1:N-2,2)...
                       + len(2:N-1)./(len(1:N-2)+len(2:N-1)).*tauelem(2:N-1,2);
    taunode(N,1:2) = tauelem(N-1,1:2);
    
    len2 = sqrt(taunode(:,1).^2+taunode(:,2).^2);
    taunode(:,1) = taunode(:,1)./len2;
    taunode(:,2) = taunode(:,2)./len2;
    
    nnode = zeros(N,2);
    nnode(:,1) = taunode(:,2);
    nnode(:,2) = -taunode(:,1);
    
    p_ss = 2*(tauelem(2:N-1,:)-tauelem(1:N-2,:))./(len(2:N-1)+len(1:N-2));
    
    kappa = zeros(N,1);
    kappa(2:N-1) = -sum(p_ss.*nnode(2:N-1,1:2),2);