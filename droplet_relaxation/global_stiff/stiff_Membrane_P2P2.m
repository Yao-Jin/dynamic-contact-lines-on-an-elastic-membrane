function [stiff_membrane] = stiff_Membrane_P2P2(geom,para,sln)

% y: P2; kappa: P2; 
% construct stiffness matrices for membrane elements
%% ----------------------------------------------------------------------------
    % read geometry data    
    N = size(geom.node,1); Ne = size(geom.edge,1);
    Nmem = para.Nmem; pmem = sln.pmem(1:Nmem,1:2); 
    mem2node = geom.mem2node; mem2edge = geom.mem2edge;
    cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2);
    
    % precondition constants on each elem
    [nelem_mem, nnode_mem, ~, ~] = getnormal_tau_mem_P1(pmem);
    pypxelem = (pmem(2:Nmem,2)-pmem(1:Nmem-1,2))./(pmem(2:Nmem,1)-pmem(1:Nmem-1,1));
    lenelem = sqrt(1+pypxelem.^2);
    pypxnode = -nnode_mem(:,1)./nnode_mem(:,2);
    lennode = sqrt(1+pypxnode.^2);
    
    kmem = sln.kmem;
    kmem2_avg = 1/2*(kmem(1:Nmem-1).^2+kmem(2:Nmem).^2);
    
%% ---------------------------------------------------------------------------
    % read model parameters
    cb = para.cb; Ca = para.Ca; gamma1 = para.gamma1; gamma2 = para.gamma2;
    
%% ---------------------------------------------------------------------------
    % initialize indices for element
%     linelemP0 = (1:Nmem-1);
%     linelemP1 = [1:Nmem-1; 2:Nmem]';
    linelemP2 = [1:Nmem-1; Nmem+1:2*Nmem-1; 2:Nmem]';
    elemP2online = [mem2node(1:end-1)' mem2edge(1:end)+N mem2node(2:end)'];
    
    i11 = zeros(9*(Nmem-1),1); j11 = zeros(9*(Nmem-1),1); s11 = zeros(9*(Nmem-1),1);
    i12 = zeros(9*(Nmem-1),1); j12 = zeros(9*(Nmem-1),1); s12 = zeros(9*(Nmem-1),1);
    i31 = zeros(9*(Nmem-1),1); j31 = zeros(9*(Nmem-1),1); s31 = zeros(9*(Nmem-1),1);
    i32 = zeros(9*(Nmem-1),1); j32 = zeros(9*(Nmem-1),1); s32 = zeros(9*(Nmem-1),1);
    i41 = zeros(9*(Nmem-1),1); j41 = zeros(9*(Nmem-1),1); s41 = zeros(9*(Nmem-1),1);
    i42 = zeros(9*(Nmem-1),1); j42 = zeros(9*(Nmem-1),1); s42 = zeros(9*(Nmem-1),1);  
    i43 = zeros(9*(Nmem-1),1); j43 = zeros(9*(Nmem-1),1); s43 = zeros(9*(Nmem-1),1);  
    i44 = zeros(9*(Nmem-1),1); j44 = zeros(9*(Nmem-1),1); s44 = zeros(9*(Nmem-1),1);
    i45 = zeros(9*(Nmem-1),1); j45 = zeros(9*(Nmem-1),1); s45 = zeros(9*(Nmem-1),1);
    i51 = zeros(9*(Nmem-1),1); j51 = zeros(9*(Nmem-1),1); s51 = zeros(9*(Nmem-1),1);
    i52 = zeros(9*(Nmem-1),1); j52 = zeros(9*(Nmem-1),1); s52 = zeros(9*(Nmem-1),1);
    
    ind = 0;
    for t = 1:Nmem-1
        
        if (t>=cl2mem1) && (t<=cl2mem2-1)
            gamma = gamma1; 
        else
            gamma = gamma2;
        end
        
        len = sqrt((pmem(t+1,1)-pmem(t,1))^2+(pmem(t+1,2)-pmem(t,2))^2);
        nelem_1 = [nnode_mem(t,1) nelem_mem(t,1) nnode_mem(t+1,1)];
        nelem_2 = [nnode_mem(t,2) nelem_mem(t,2) nnode_mem(t+1,2)];

        len_elem = [lennode(t) lenelem(t) lennode(t+1)];

        Ls11 = 3*cb/(2*Ca)*kmem2_avg(t)/len*getStiffonRef1DP2der; Ls11 = Ls11.*(len_elem.*nelem_1)';
        Ls12 = 3*cb/(2*Ca)*kmem2_avg(t)/len*getStiffonRef1DP2der; Ls12 = Ls12.*(len_elem.*nelem_2)';     
        
        Ls31 = cb/Ca/len*getStiffonRef1DP2der; Ls31 = Ls31.*(len_elem.*nelem_1)'.*(1./len_elem);
        Ls32 = cb/Ca/len*getStiffonRef1DP2der; Ls32 = Ls32.*(len_elem.*nelem_2)'.*(1./len_elem);       
        
        Ls41 = len*getStiffonRef1DP2P2;
        Ls42 = len*getStiffonRef1DP2P2; Ls42 = Ls42.*(len_elem.*nelem_1);
        Ls43 = len*getStiffonRef1DP2P2; Ls43 = Ls43.*(len_elem.*nelem_2);
        Ls44 = len*getStiffonRef1DP2P2; Ls44 = Ls44./len_elem;
        Ls45 = 1/len*getStiffonRef1DP2der; 
        
        Ls51 = gamma/Ca/len*getStiffonRef1DP2der; Ls51 = Ls51.*(len_elem.*nelem_1)';
        Ls52 = gamma/Ca/len*getStiffonRef1DP2der; Ls52 = Ls52.*(len_elem.*nelem_2)';
        
        for ti = 1:3
            for tj = 1:3
                ind = ind + 1;
                i11(ind) = elemP2online(t,ti); j11(ind) = linelemP2(t,tj); s11(ind) = Ls11(ti,tj);
                i12(ind) = elemP2online(t,ti); j12(ind) = linelemP2(t,tj); s12(ind) = Ls12(ti,tj);
                i31(ind) = elemP2online(t,ti); j31(ind) = linelemP2(t,tj); s31(ind) = Ls31(ti,tj);
                i32(ind) = elemP2online(t,ti); j32(ind) = linelemP2(t,tj); s32(ind) = Ls32(ti,tj);               
                i41(ind) = linelemP2(t,ti);    j41(ind) = linelemP2(t,tj); s41(ind) = Ls41(ti,tj);      
                i42(ind) = linelemP2(t,tj);    j42(ind) = elemP2online(t,ti); s42(ind) = Ls42(tj,ti);
                i43(ind) = linelemP2(t,tj);    j43(ind) = elemP2online(t,ti); s43(ind) = Ls43(tj,ti);
                i44(ind) = linelemP2(t,ti);    j44(ind) = linelemP2(t,tj); s44(ind) = Ls44(ti,tj);
                i45(ind) = linelemP2(t,ti);    j45(ind) = linelemP2(t,tj); s45(ind) = Ls45(ti,tj);
                
                i51(ind) = elemP2online(t,ti); j51(ind) = linelemP2(t,tj); s51(ind) = Ls51(ti,tj);
                i52(ind) = elemP2online(t,ti); j52(ind) = linelemP2(t,tj); s52(ind) = Ls52(ti,tj);
            end
        end
        
    end
    
    stiff_membrane = [];
    stiff_membrane.M11 = sparse(i11,j11,s11,N+Ne,2*Nmem-1);
    stiff_membrane.M12 = sparse(i12,j12,s12,N+Ne,2*Nmem-1);
    
    stiff_membrane.M31 = sparse(i31,j31,s31,N+Ne,2*Nmem-1);
    stiff_membrane.M32 = sparse(i32,j32,s32,N+Ne,2*Nmem-1);  
    
    stiff_membrane.M41 = sparse(i41,j41,s41,2*Nmem-1,2*Nmem-1);
    stiff_membrane.M42 = sparse(i42,j42,s42,2*Nmem-1,N+Ne);
    stiff_membrane.M43 = sparse(i43,j43,s43,2*Nmem-1,N+Ne);
    stiff_membrane.M44 = sparse(i44,j44,s44,2*Nmem-1,2*Nmem-1);
    stiff_membrane.M45 = sparse(i45,j45,s45,2*Nmem-1,2*Nmem-1);

    stiff_membrane.M51 = sparse(i51,j51,s51,N+Ne,2*Nmem-1);
    stiff_membrane.M52 = sparse(i52,j52,s52,N+Ne,2*Nmem-1);
    
end