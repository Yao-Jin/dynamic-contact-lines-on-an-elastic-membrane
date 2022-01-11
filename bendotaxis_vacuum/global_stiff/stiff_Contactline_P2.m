function [stiff_cl] = stiff_Contactline_P2(geom,para,sln)

% construct stiffness matrices or load vectors for elements on contact lines
%% read geometry parameters
    N = size(geom.node,1); Ne = size(geom.edge,1); Nmem = para.Nmem;
    cl2node1 = geom.cl2node(1); cl2node2 = geom.cl2node(2);
    cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2); Nm = cl2mem2-cl2mem1+1;
    pmem = sln.pmem; x = pmem(:,1); y = pmem(:,2);
    [~, nnode_mem, ~, taunode_mem] = getnormal_tau_mem_P1(pmem);
    pypxelem = (pmem(2:Nmem,2)-pmem(1:Nmem-1,2))./(pmem(2:Nmem,1)-pmem(1:Nmem-1,1)); pypxnode = -nnode_mem(:,1)./nnode_mem(:,2);
    lenelem = sqrt(1+pypxelem.^2);
    lennode = -1./nnode_mem(:,2);
    
%% read model parameters
    Ca = para.Ca; gamma1 = para.gamma1; gamma2 = para.gamma2; cb = para.cb; mu0 = para.mu0; %ls = para.ls;
%% --------------------------------------------------------------------------        
    stiff_cl = [];
    
    val1 = -(gamma1-gamma2)/Ca*taunode_mem(cl2mem1,1);
    val2 = (gamma1-gamma2)/Ca*taunode_mem(cl2mem2,1);
    val3 = -(gamma1-gamma2)/Ca*taunode_mem(cl2mem1,2);
    val4 = (gamma1-gamma2)/Ca*taunode_mem(cl2mem2,2);
    stiff_cl.cl1 = sparse([cl2node1 cl2node2],[1 1],[val1 val2],N+Ne,1);
    stiff_cl.cl2 = sparse([cl2node1 cl2node2],[1 1],[val3 val4],N+Ne,1);
    
    stiff_cl.cl3 = sparse(cl2node1,Nm+1,taunode_mem(cl2mem1,1)/Ca,N+Ne,2*Nm);
    stiff_cl.cl4 = sparse(cl2node1,Nm+1,taunode_mem(cl2mem1,2)/Ca,N+Ne,2*Nm);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
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

    stiff_cl.CL1 = sparse([cl2node1 cl2node2],[1 1],(gamma1-gamma2)/Ca*[-loclen1 loclen2],N+Ne,1);    
    
    con1 = -cb/Ca/loclen1; con2 = cb/Ca/loclen2;

%----------------------------------------------------------------------------
    cl5 = sparse(N+Ne,2*Nmem-1); cl6 = sparse(N+Ne,2*Nmem-1);
    
    dxL = x(cl2mem1)-x(cl2mem1-1);
    cl5(cl2node1,cl2mem1) = 3/dxL*con1*nnode_mem(cl2mem1,1);
    cl5(cl2node1,Nmem+cl2mem1-1) = -4/dxL*con1*nnode_mem(cl2mem1,1);
    cl5(cl2node1,cl2mem1-1) = 1/dxL*con1*nnode_mem(cl2mem1,1);
    
    cl6(cl2node1,cl2mem1) = 3/dxL*con1*nnode_mem(cl2mem1,2);
    cl6(cl2node1,Nmem+cl2mem1-1) = -4/dxL*con1*nnode_mem(cl2mem1,2);
    cl6(cl2node1,cl2mem1-1) = 1/dxL*con1*nnode_mem(cl2mem1,2);    
    
    dxR = x(cl2mem2+1)-x(cl2mem2);
    cl5(cl2node2,cl2mem2+1) = -1/dxR*con2*nnode_mem(cl2mem2,1);
    cl5(cl2node2,Nmem+cl2mem2) = 4/dxR*con2*nnode_mem(cl2mem2,1);
    cl5(cl2node2,cl2mem2) = -3/dxR*con2*nnode_mem(cl2mem2,1);
    
    cl6(cl2node2,cl2mem2+1) = -1/dxR*con2*nnode_mem(cl2mem2,2);
    cl6(cl2node2,Nmem+cl2mem2) = 4/dxR*con2*nnode_mem(cl2mem2,2);
    cl6(cl2node2,cl2mem2) = -3/dxR*con2*nnode_mem(cl2mem2,2);        
%----------------------------------------------------------------------------
    stiff_cl.cl5 = cl5;
    stiff_cl.cl6 = cl6;   
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    stiff_cl.cl7 = sparse([1 Nm 1],[1 Nm Nm+1],[1/mu0 1/mu0 -1/mu0],2*Nm,2*Nm);
    
    kmem = sln.kmem;
    linelemP2 = [1:Nmem-1; Nmem+1:2*Nmem-1; 2:Nmem]';
    Nl = cl2mem1;
    i = zeros(3*(Nl-1),1); j = zeros(3*(Nl-1),1); s = zeros(3*(Nl-1),1);
    ind = 0;
    for t = 1:cl2mem1-1
        len = sqrt((pmem(t+1,1)-pmem(t,1))^2+(pmem(t+1,2)-pmem(t,2))^2);
        len_elem = [lennode(t) lenelem(t) lennode(t+1)];
        Ls = len*[kmem(t) kmem(Nmem+t) kmem(t+1)]./len_elem.*[1/6 2/3 1/6];
        for k = 1:3
            ind = ind + 1;
            i(ind) = Nm+1; j(ind) = linelemP2(t,k); s(ind) = Ls(k);
        end
    end
    stiff_cl.cl8 = sparse(i,j,s,2*Nm,2*Nmem-1);
    
    stiff_cl.cl9 = sparse(Nm+1,cl2node1,taunode_mem(cl2mem1,1),2*Nm,N+Ne);
    stiff_cl.cl10 = sparse(Nm+1,cl2node1,taunode_mem(cl2mem1,2),2*Nm,N+Ne);

    stiff_cl.cl11 = sparse([Nm+1 Nm+1],[Nm+1 1],[1/mu0/Ca -1/mu0/Ca],2*Nm,2*Nm);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    val1 = -(gamma1-gamma2)/Ca*nnode_mem(cl2mem1,1)*pypxnode(cl2mem1);
    val2 = (gamma1-gamma2)/Ca*nnode_mem(cl2mem2,1)*pypxnode(cl2mem2);
    val3 = -(gamma1-gamma2)/Ca*nnode_mem(cl2mem1,2)*pypxnode(cl2mem1);
    val4 = (gamma1-gamma2)/Ca*nnode_mem(cl2mem2,2)*pypxnode(cl2mem2);
    stiff_cl.cl12 = sparse([cl2node1 cl2node2],[1 1],[val1 val2],N+Ne,1);
    stiff_cl.cl13 = sparse([cl2node1 cl2node2],[1 1],[val3 val4],N+Ne,1);
end
    