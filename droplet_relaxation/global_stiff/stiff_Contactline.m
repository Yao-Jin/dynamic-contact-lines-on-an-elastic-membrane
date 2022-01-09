function [stiff_cl] = stiff_Contactline(geom,para,sln)

% construct stiffness matrices or load vectors for elements on contact lines
%% read geometry parameters
    N = size(geom.node,1); Ne = size(geom.edge,1);
    cl2node1 = geom.cl2node(1); cl2node2 = geom.cl2node(2);
    cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2);
    Nmem = para.Nmem;
%     pmem = sln.pmem; pmem_nodal = pmem(1:Nmem,1:2);
    pmem_nodal = sln.pmem(1:Nmem,1:2); x = pmem_nodal(:,1); y = pmem_nodal(:,2);
    [~, ~, ~, taunode_mem] = getnormal_tau_mem_P1(pmem_nodal);
    
%% read model parameters
    Ca = para.Ca; gamma1 = para.gamma1; gamma2 = para.gamma2;
    mu0 = para.mu0; ls = para.ls;
%% --------------------------------------------------------------------------        
    stiff_cl = [];
    
    pypxcl1 = (x(cl2mem1+1)-x(cl2mem1))/(x(cl2mem1+1)-x(cl2mem1-1))*...
                    (y(cl2mem1)-y(cl2mem1-1))/(x(cl2mem1)-x(cl2mem1-1))...
            + (x(cl2mem1)-x(cl2mem1-1))/(x(cl2mem1+1)-x(cl2mem1-1))*...
                    (y(cl2mem1+1)-y(cl2mem1))/(x(cl2mem1+1)-x(cl2mem1));
    loclen1 = sqrt(1+pypxcl1^2);            
    pypxcl2 = (x(cl2mem2+1)-x(cl2mem2))/(x(cl2mem2+1)-x(cl2mem2-1))*...
                    (y(cl2mem2)-y(cl2mem2-1))/(x(cl2mem2)-x(cl2mem2-1))...
            + (x(cl2mem2)-x(cl2mem2-1))/(x(cl2mem2+1)-x(cl2mem2-1))*...
                    (y(cl2mem2+1)-y(cl2mem2))/(x(cl2mem2+1)-x(cl2mem2));
    loclen2 = sqrt(1+pypxcl2^2);      
    
    val1 = -(gamma1-gamma2)/Ca*loclen1;
    val2 = (gamma1-gamma2)/Ca*loclen2;

    stiff_cl.cl1 = sparse([cl2node1 cl2node2],[1 1],[val1 val2],N+Ne,1);

%% --------------------------------------------------------------------------  
    i = [cl2mem1 cl2mem1 Nmem+1 Nmem+1 cl2mem2 cl2mem2 Nmem+2 Nmem+2];
    j = [cl2mem1 Nmem+1 cl2mem1 Nmem+1 cl2mem2 Nmem+2 cl2mem2 Nmem+2];
    s = 1/mu0/ls*[1 -1 -1 1 1 -1 -1 1];
    stiff_cl.cl4 = sparse(i,j,s,2*Nmem+1,2*Nmem+1);  
    
    val1 = -(gamma1-gamma2)/Ca*taunode_mem(cl2mem1,1);
    val2 = (gamma1-gamma2)/Ca*taunode_mem(cl2mem2,1);
    val3 = -(gamma1-gamma2)/Ca*taunode_mem(cl2mem1,2);
    val4 = (gamma1-gamma2)/Ca*taunode_mem(cl2mem2,2);
    stiff_cl.cl5 = sparse([cl2node1 cl2node2],[1 1],[val1 val2],N+Ne,1);
    stiff_cl.cl6 = sparse([cl2node1 cl2node2],[1 1],[val3 val4],N+Ne,1);    
        
end
    