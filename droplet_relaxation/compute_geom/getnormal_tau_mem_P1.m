function [nelem, nnode, tauelem, taunode] = getnormal_tau_mem_P1(pnodal)

    N = size(pnodal,1); x = pnodal(:,1); y = pnodal(:,2);
    
    Dy_elem = (y(2:end)-y(1:end-1))./(x(2:end)-x(1:end-1));
    len = sqrt(1+Dy_elem.^2);
    tauelem = [1./len Dy_elem./len];
    nelem = [tauelem(:,2) -tauelem(:,1)];
    
    a = [(x(3)-x(1))/(x(3)-x(2)); (x(3:end)-x(2:end-1))./(x(3:end)-x(1:end-2)); (x(N)-x(N-2))/(x(N-1)-x(N-2))];
    b = [(x(1)-x(2))/(x(3)-x(2)); (x(2:end-1)-x(1:end-2))./(x(3:end)-x(1:end-2)); (x(N-1)-x(N))/(x(N-1)-x(N-2))];
 
    Dy_node = a.*[Dy_elem(1); Dy_elem(1:end-1); Dy_elem(N-1)]...
            + b.*[(y(3)-y(1))/(x(3)-x(1)); Dy_elem(2:end); (y(N)-y(N-2))/(x(N)-x(N-2))];
    
    len = sqrt(1+Dy_node.^2);
    taunode = [1./len Dy_node./len];
    
    nnode = [taunode(:,2) -taunode(:,1)];

end