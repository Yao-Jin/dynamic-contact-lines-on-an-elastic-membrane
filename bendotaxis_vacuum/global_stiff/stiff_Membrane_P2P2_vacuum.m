function [stiff_membrane] = stiff_Membrane_P2P2_vacuum(geom,para,sln)

% y: P1; kappa: P1; 
% construct stiffness matrices for membrane elements
%% ----------------------------------------------------------------------------
    % read geometry data    
    N = size(geom.node,1); Ne = size(geom.edge,1); Nmem = para.Nmem; 
    pmem = sln.pmem(1:Nmem,1:2); mem2node = geom.mem2node; mem2edge = geom.mem2edge;
    cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2); 
    Nl = cl2mem1; Nm = cl2mem2-cl2mem1+1; Nr = Nmem-cl2mem2+1;
    
    % precondition constants on each elem
    [nelem_mem, nnode_mem, ~, ~] = getnormal_tau_mem_P1(pmem);
    pypxelem = (pmem(2:Nmem,2)-pmem(1:Nmem-1,2))./(pmem(2:Nmem,1)-pmem(1:Nmem-1,1));
    lenelem = sqrt(1+pypxelem.^2);
    lennode = -1./nnode_mem(:,2);
    
    kmem = sln.kmem;
    kmem2_avg = 1/2*(kmem(1:Nmem-1).^2+kmem(2:Nmem).^2);
    
%% ---------------------------------------------------------------------------
    % read model parameters
    cb = para.cb; Ca = para.Ca; gamma1 = para.gamma1; gamma2 = para.gamma2;
    
%% ---------------------------------------------------------------------------
    % initialize indices for element
    linelemP2 = [1:Nmem-1; Nmem+1:2*Nmem-1; 2:Nmem]';
    elemP2online = [mem2node(1:end-1) mem2edge(1:end)+N mem2node(2:end)];
    
%% -----------------------------------------------------------------------------------
    % membrane under the droplet
    % sigma_1
    i11 = zeros(9*(Nm-1),1); j11 = zeros(9*(Nm-1),1); s11 = zeros(9*(Nm-1),1);
    i12 = zeros(9*(Nm-1),1); j12 = zeros(9*(Nm-1),1); s12 = zeros(9*(Nm-1),1);
    i13 = zeros(9*(Nm-1),1); j13 = zeros(9*(Nm-1),1); s13 = zeros(9*(Nm-1),1);
    i14 = zeros(9*(Nm-1),1); j14 = zeros(9*(Nm-1),1); s14 = zeros(9*(Nm-1),1);    
    i15 = zeros(9*(Nm-1),1); j15 = zeros(9*(Nm-1),1); s15 = zeros(9*(Nm-1),1);
    i16 = zeros(9*(Nm-1),1); j16 = zeros(9*(Nm-1),1); s16 = zeros(9*(Nm-1),1);
    i17 = zeros(9*(Nm-1),1); j17 = zeros(9*(Nm-1),1); s17 = zeros(9*(Nm-1),1);
    i18 = zeros(9*(Nm-1),1); j18 = zeros(9*(Nm-1),1); s18 = zeros(9*(Nm-1),1);
    i19 = zeros(9*(Nm-1),1); j19 = zeros(9*(Nm-1),1); s19 = zeros(9*(Nm-1),1);
    
    ind = 0;
    for t = cl2mem1:cl2mem2-1
        
        len = sqrt((pmem(t+1,1)-pmem(t,1))^2+(pmem(t+1,2)-pmem(t,2))^2);
        nelem_1 = [nnode_mem(t,1) nelem_mem(t,1) nnode_mem(t+1,1)];
        nelem_2 = [nnode_mem(t,2) nelem_mem(t,2) nnode_mem(t+1,2)];
        len_elem = [lennode(t) lenelem(t) lennode(t+1)];
        
        Ls11 = gamma1/Ca/len*getStiffonRef1DP2der; Ls11 = Ls11.*(len_elem.*nelem_1)';
        Ls12 = gamma1/Ca/len*getStiffonRef1DP2der; Ls12 = Ls12.*(len_elem.*nelem_2)';
        
        Ls13 = cb/Ca/len*getStiffonRef1DP2der; Ls13 = Ls13.*(len_elem.*nelem_1)'.*(1./len_elem);
        Ls14 = cb/Ca/len*getStiffonRef1DP2der; Ls14 = Ls14.*(len_elem.*nelem_2)'.*(1./len_elem);
        
        Ls15 = 3*cb/(2*Ca)*kmem2_avg(t)/len*getStiffonRef1DP2der; Ls15 = Ls15.*(len_elem.*nelem_1)';
        Ls16 = 3*cb/(2*Ca)*kmem2_avg(t)/len*getStiffonRef1DP2der; Ls16 = Ls16.*(len_elem.*nelem_2)';      
        
        Ls17 = len*getStiffonRef1DP2P2;
        Ls18 = len*getStiffonRef1DP2P2; Ls18 = Ls18.*(len_elem.*nelem_1);
        Ls19 = len*getStiffonRef1DP2P2; Ls19 = Ls19.*(len_elem.*nelem_2);
                
        tt = t-cl2mem1+1;
        for ti = 1:3
            for tj = 1:3
                ind = ind + 1;
                i11(ind) = elemP2online(tt,ti); j11(ind) = linelemP2(t,tj); s11(ind) = Ls11(ti,tj);
                i12(ind) = elemP2online(tt,ti); j12(ind) = linelemP2(t,tj); s12(ind) = Ls12(ti,tj);
                i13(ind) = elemP2online(tt,ti); j13(ind) = linelemP2(t,tj); s13(ind) = Ls13(ti,tj);
                i14(ind) = elemP2online(tt,ti); j14(ind) = linelemP2(t,tj); s14(ind) = Ls14(ti,tj);
                i15(ind) = elemP2online(tt,ti); j15(ind) = linelemP2(t,tj); s15(ind) = Ls15(ti,tj);
                i16(ind) = elemP2online(tt,ti); j16(ind) = linelemP2(t,tj); s16(ind) = Ls16(ti,tj);
                
                i17(ind) = linelemP2(t,ti); j17(ind) = linelemP2(t,tj); s17(ind) = Ls17(ti,tj);
                i18(ind) = linelemP2(t,ti); j18(ind) = elemP2online(tt,tj); s18(ind) = Ls18(ti,tj);
                i19(ind) = linelemP2(t,ti); j19(ind) = elemP2online(tt,tj); s19(ind) = Ls19(ti,tj);                
            end
        end       
    end
    
% --------------------------------------------------------------------------------------------------    
    t = cl2mem1-1;    len = sqrt((pmem(t+1,1)-pmem(t,1))^2+(pmem(t+1,2)-pmem(t,2))^2);
    nelem_1 = [nnode_mem(t,1) nelem_mem(t,1) nnode_mem(t+1,1)];
    nelem_2 = [nnode_mem(t,2) nelem_mem(t,2) nnode_mem(t+1,2)];
    len_elem = [lennode(t) lenelem(t) lennode(t+1)];
        
    Ls11 = gamma2/Ca/len*getStiffonRef1DP2der; Ls11 = Ls11.*(len_elem.*nelem_1)';
    Ls12 = gamma2/Ca/len*getStiffonRef1DP2der; Ls12 = Ls12.*(len_elem.*nelem_2)';
        
    Ls13 = cb/Ca/len*getStiffonRef1DP2der; Ls13 = Ls13.*(len_elem.*nelem_1)'.*(1./len_elem);
    Ls14 = cb/Ca/len*getStiffonRef1DP2der; Ls14 = Ls14.*(len_elem.*nelem_2)'.*(1./len_elem);
        
    Ls15 = 3*cb/(2*Ca)*kmem2_avg(t)/len*getStiffonRef1DP2der; Ls15 = Ls15.*(len_elem.*nelem_1)';
    Ls16 = 3*cb/(2*Ca)*kmem2_avg(t)/len*getStiffonRef1DP2der; Ls16 = Ls16.*(len_elem.*nelem_2)'; 
    
    for tj = 1:3
        ind = ind + 1;
        i11(ind) = mem2node(1); j11(ind) = linelemP2(t,tj); s11(ind) = Ls11(3,tj);
        i12(ind) = mem2node(1); j12(ind) = linelemP2(t,tj); s12(ind) = Ls12(3,tj);
        i13(ind) = mem2node(1); j13(ind) = linelemP2(t,tj); s13(ind) = Ls13(3,tj);
        i14(ind) = mem2node(1); j14(ind) = linelemP2(t,tj); s14(ind) = Ls14(3,tj);
        i15(ind) = mem2node(1); j15(ind) = linelemP2(t,tj); s15(ind) = Ls15(3,tj);
        i16(ind) = mem2node(1); j16(ind) = linelemP2(t,tj); s16(ind) = Ls16(3,tj);               
    end       
    
% --------------------------------------------------------------------------------------------------    
    t = cl2mem2;    len = sqrt((pmem(t+1,1)-pmem(t,1))^2+(pmem(t+1,2)-pmem(t,2))^2);
    nelem_1 = [nnode_mem(t,1) nelem_mem(t,1) nnode_mem(t+1,1)];
    nelem_2 = [nnode_mem(t,2) nelem_mem(t,2) nnode_mem(t+1,2)];
    len_elem = [lennode(t) lenelem(t) lennode(t+1)];
        
    Ls11 = gamma2/Ca/len*getStiffonRef1DP2der; Ls11 = Ls11.*(len_elem.*nelem_1)';
    Ls12 = gamma2/Ca/len*getStiffonRef1DP2der; Ls12 = Ls12.*(len_elem.*nelem_2)';
        
    Ls13 = cb/Ca/len*getStiffonRef1DP2der; Ls13 = Ls13.*(len_elem.*nelem_1)'.*(1./len_elem);
    Ls14 = cb/Ca/len*getStiffonRef1DP2der; Ls14 = Ls14.*(len_elem.*nelem_2)'.*(1./len_elem);
        
    Ls15 = 3*cb/(2*Ca)*kmem2_avg(t)/len*getStiffonRef1DP2der; Ls15 = Ls15.*(len_elem.*nelem_1)';
    Ls16 = 3*cb/(2*Ca)*kmem2_avg(t)/len*getStiffonRef1DP2der; Ls16 = Ls16.*(len_elem.*nelem_2)'; 
    
    for tj = 1:3
        ind = ind + 1;
        i11(ind) = mem2node(end); j11(ind) = linelemP2(t,tj); s11(ind) = Ls11(1,tj);
        i12(ind) = mem2node(end); j12(ind) = linelemP2(t,tj); s12(ind) = Ls12(1,tj);
        i13(ind) = mem2node(end); j13(ind) = linelemP2(t,tj); s13(ind) = Ls13(1,tj);
        i14(ind) = mem2node(end); j14(ind) = linelemP2(t,tj); s14(ind) = Ls14(1,tj);
        i15(ind) = mem2node(end); j15(ind) = linelemP2(t,tj); s15(ind) = Ls15(1,tj);
        i16(ind) = mem2node(end); j16(ind) = linelemP2(t,tj); s16(ind) = Ls16(1,tj);               
    end 
% --------------------------------------------------------------------------------------------------    
    
    stiff_membrane = [];
    stiff_membrane.M11 = sparse(i11,j11,s11,N+Ne,2*Nmem-1);
    stiff_membrane.M12 = sparse(i12,j12,s12,N+Ne,2*Nmem-1);
    stiff_membrane.M13 = sparse(i13,j13,s13,N+Ne,2*Nmem-1);
    stiff_membrane.M14 = sparse(i14,j14,s14,N+Ne,2*Nmem-1);
    stiff_membrane.M15 = sparse(i15,j15,s15,N+Ne,2*Nmem-1);
    stiff_membrane.M16 = sparse(i16,j16,s16,N+Ne,2*Nmem-1);
    stiff_membrane.M17 = sparse(i17,j17,s17,2*Nmem-1,2*Nmem-1);
    stiff_membrane.M18 = sparse(i18,j18,s18,2*Nmem-1,N+Ne);
    stiff_membrane.M19 = sparse(i19,j19,s19,2*Nmem-1,N+Ne);    

%% -----------------------------------------------------------------------------------
    % membrane sigma_2l
    i21 = zeros(9*(Nl-1),1); j21 = zeros(9*(Nl-1),1); s21 = zeros(9*(Nl-1),1);
    i22 = zeros(9*(Nl-1),1); j22 = zeros(9*(Nl-1),1); s22 = zeros(9*(Nl-1),1);    
    i23 = zeros(9*(Nl-1),1); j23 = zeros(9*(Nl-1),1); s23 = zeros(9*(Nl-1),1); 
    i24 = zeros(9*(Nl-1),1); j24 = zeros(9*(Nl-1),1); s24 = zeros(9*(Nl-1),1);
    ind = 0;
    for t = 1:cl2mem1-1
        len = sqrt((pmem(t+1,1)-pmem(t,1))^2+(pmem(t+1,2)-pmem(t,2))^2);
        len_elem = [lennode(t) lenelem(t) lennode(t+1)];
        
        Ls21 = gamma2/len*getStiffonRef1DP2der; Ls21 = Ls21.*(len_elem)';
        Ls22 = cb/len*getStiffonRef1DP2der; Ls22 = Ls22.*(len_elem)'.*(1./len_elem);
        Ls23 = 3*cb/2*kmem2_avg(t)/len*getStiffonRef1DP2der; Ls23 = Ls23.*(len_elem)';
        Ls24 = len*getStiffonRef1DP2P2;
        for ti = 1:3
            for tj = 1:3
                ind = ind + 1;
                i21(ind) = linelemP2(t,ti); j21(ind) = linelemP2(t,tj); s21(ind) = Ls21(ti,tj);    
                i22(ind) = linelemP2(t,ti); j22(ind) = linelemP2(t,tj); s22(ind) = Ls22(ti,tj); 
                i23(ind) = linelemP2(t,ti); j23(ind) = linelemP2(t,tj); s23(ind) = Ls23(ti,tj);
                i24(ind) = linelemP2(t,ti); j24(ind) = linelemP2(t,tj); s24(ind) = Ls24(ti,tj);
            end
        end        
    end
    
    M21 = sparse(i21,j21,s21,2*Nmem-1,2*Nmem-1);
    stiff_membrane.M21 = [M21(1:cl2mem1-1,:); sparse(1,2*Nmem-1); M21(cl2mem1+1:2*Nmem-1,:)];
    M22 = sparse(i22,j22,s22,2*Nmem-1,2*Nmem-1);    
    stiff_membrane.M22 = [M22(1:cl2mem1-1,:); sparse(1,2*Nmem-1); M22(cl2mem1+1:2*Nmem-1,:)];
    M23 = sparse(i23,j23,s23,2*Nmem-1,2*Nmem-1);    
    stiff_membrane.M23 = [M23(1:cl2mem1-1,:); sparse(1,2*Nmem-1); M23(cl2mem1+1:2*Nmem-1,:)]; 
    M24 = sparse(i24,j24,s24,2*Nmem-1,2*Nmem-1);    
    stiff_membrane.M24 = [M24(1:cl2mem1-1,:); sparse(1,2*Nmem-1); M24(cl2mem1+1:2*Nmem-1,:)]; 

%% -----------------------------------------------------------------------------------
    % membrane sigma_2r
    i31 = zeros(9*(Nr-1),1); j31 = zeros(9*(Nr-1),1); s31 = zeros(9*(Nr-1),1);
    i32 = zeros(9*(Nr-1),1); j32 = zeros(9*(Nr-1),1); s32 = zeros(9*(Nr-1),1);   
    i33 = zeros(9*(Nr-1),1); j33 = zeros(9*(Nr-1),1); s33 = zeros(9*(Nr-1),1); 
    ind = 0;
    for t = cl2mem2:Nmem-1
        len = sqrt((pmem(t+1,1)-pmem(t,1))^2+(pmem(t+1,2)-pmem(t,2))^2);
        len_elem = [lennode(t) lenelem(t) lennode(t+1)];
        
        Ls31 = gamma2/len*getStiffonRef1DP2der; Ls31 = Ls31.*(len_elem)';
        Ls32 = cb/len*getStiffonRef1DP2der; Ls32 = Ls32.*(len_elem)'.*(1./len_elem);
        Ls33 = 3*cb/2*kmem2_avg(t)/len*getStiffonRef1DP2der; Ls33 = Ls33.*(len_elem)';

        for ti = 1:3
            for tj = 1:3
                ind = ind + 1;
                i31(ind) = linelemP2(t,ti); j31(ind) = linelemP2(t,tj); s31(ind) = Ls31(ti,tj);    
                i32(ind) = linelemP2(t,ti); j32(ind) = linelemP2(t,tj); s32(ind) = Ls32(ti,tj); 
                i33(ind) = linelemP2(t,ti); j33(ind) = linelemP2(t,tj); s33(ind) = Ls33(ti,tj); 
            end
        end        
    end
    
    M31 = sparse(i31,j31,s31,2*Nmem-1,2*Nmem-1);
    stiff_membrane.M31 = [M31(1:cl2mem2-1,:); sparse(1,2*Nmem-1); M31(cl2mem2+1:2*Nmem-1,:)];
    M32 = sparse(i32,j32,s32,2*Nmem-1,2*Nmem-1);    
    stiff_membrane.M32 = [M32(1:cl2mem2-1,:); sparse(1,2*Nmem-1); M32(cl2mem2+1:2*Nmem-1,:)];
    M33 = sparse(i33,j33,s33,2*Nmem-1,2*Nmem-1);    
    stiff_membrane.M33 = [M33(1:cl2mem2-1,:); sparse(1,2*Nmem-1); M33(cl2mem2+1:2*Nmem-1,:)];     
%% -----------------------------------------------------------------------------------
    % whole membrane
    i41 = zeros(9*(Nmem-1),1); j41 = zeros(9*(Nmem-1),1); s41 = zeros(9*(Nmem-1),1);
    i42 = zeros(9*(Nmem-1),1); j42 = zeros(9*(Nmem-1),1); s42 = zeros(9*(Nmem-1),1);    
    ind = 0;
    for t = 1:Nmem-1
        len = sqrt((pmem(t+1,1)-pmem(t,1))^2+(pmem(t+1,2)-pmem(t,2))^2);
        len_elem = [lennode(t) lenelem(t) lennode(t+1)];
        
        Ls41 = len*getStiffonRef1DP2P2; Ls41 = Ls41./len_elem;
        Ls42 = 1/len*getStiffonRef1DP2der; 
        
        for ti = 1:3
            for tj = 1:3
                ind = ind + 1;
                i41(ind) = linelemP2(t,ti); j41(ind) = linelemP2(t,tj); s41(ind) = Ls41(ti,tj);    
                i42(ind) = linelemP2(t,ti); j42(ind) = linelemP2(t,tj); s42(ind) = Ls42(ti,tj); 
            end
        end        
    end
    
    stiff_membrane.M41 = sparse(i41,j41,s41,2*Nmem-1,2*Nmem-1);
    stiff_membrane.M42 = sparse(i42,j42,s42,2*Nmem-1,2*Nmem-1);    
end