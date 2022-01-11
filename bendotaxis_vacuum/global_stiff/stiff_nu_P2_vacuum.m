function [stiff_nu] = stiff_nu_P2_vacuum(geom,para,sln)

% nu: P1;
% construct stiffness matrices for membrane elements
%% ----------------------------------------------------------------------------
    % read geometry data
    N = size(geom.node,1); Ne = size(geom.edge,1); Nmem = para.Nmem;
    pmem = sln.pmem; mem2node = geom.mem2node; mem2edge = geom.mem2edge;
    cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2);
    Nl = cl2mem1; Nm = cl2mem2-cl2mem1+1;
    
    % precondition constants on each elem
    [nelem_mem, nnode_mem, tauelem_mem, ~] = getnormal_tau_mem_P1(pmem);
    pypxelem = (pmem(2:Nmem,2)-pmem(1:Nmem-1,2))./(pmem(2:Nmem,1)-pmem(1:Nmem-1,1));
    lenelem = sqrt(1+pypxelem.^2);
    
    kmem = sln.kmem;
    kmem_avg = 1/2*(kmem(1:Nmem-1)+kmem(2:Nmem));
%% ---------------------------------------------------------------------------
    % read model parameters
    ls = para.ls; mu = para.mu1; Ca = para.Ca;

%% ---------------------------------------------------------------------------
    % initialize indices for element
%     linelemP1 = [1:Nmem-1; 2:Nmem]';
    linelemP2 = [1:Nm-1; Nm+2:2*Nm; 2:Nm]';
    elemP2online = [mem2node(1:end-1) mem2edge(1:end)+N mem2node(2:end)];

    i11 = zeros(9*(Nm-1),1); j11 = zeros(9*(Nm-1),1); s11 = zeros(9*(Nm-1),1);
    i12 = zeros(9*(Nm-1),1); j12 = zeros(9*(Nm-1),1); s12 = zeros(9*(Nm-1),1);
    i13 = zeros(9*(Nm-1),1); j13 = zeros(9*(Nm-1),1); s13 = zeros(9*(Nm-1),1);
    i14 = zeros(9*(Nm-1),1); j14 = zeros(9*(Nm-1),1); s14 = zeros(9*(Nm-1),1);
    i15 = zeros(9*(Nm-1),1); j15 = zeros(9*(Nm-1),1); s15 = zeros(9*(Nm-1),1);
    i16 = zeros(9*(Nm-1),1); j16 = zeros(9*(Nm-1),1); s16 = zeros(9*(Nm-1),1);
    i17 = zeros(9*(Nm-1),1); j17 = zeros(9*(Nm-1),1); s17 = zeros(9*(Nm-1),1);
    
    ind = 0;
    for t = cl2mem1:cl2mem2-1
        Ls11 = 1/Ca*tauelem_mem(t,1)*getStiffonRef1DP2P2der;
        Ls12 = 1/Ca*tauelem_mem(t,2)*getStiffonRef1DP2P2der;
        Ls13 = ls/mu/((pmem(t+1,1)-pmem(t,1))*lenelem(t))*getStiffonRef1DP2der;
        Ls14 = Ca*tauelem_mem(t,1)*getStiffonRef1DP2P2der;
        Ls15 = Ca*tauelem_mem(t,2)*getStiffonRef1DP2P2der;
        
        Ls16 = kmem_avg(t)*nelem_mem(t,1)*((pmem(t+1,1)-pmem(t,1))*lenelem(t))*getStiffonRef1DP2P2;
        Ls17 = kmem_avg(t)*nelem_mem(t,2)*((pmem(t+1,1)-pmem(t,1))*lenelem(t))*getStiffonRef1DP2P2;
        tt = t-cl2mem1+1;
        for ti = 1:3
            for tj = 1:3
                ind = ind + 1;
                i11(ind) = elemP2online(tt,tj); j11(ind) = linelemP2(tt,ti); s11(ind) = Ls11(ti,tj);
                i12(ind) = elemP2online(tt,tj); j12(ind) = linelemP2(tt,ti); s12(ind) = Ls12(ti,tj);
                i13(ind) = linelemP2(tt,ti); j13(ind) = linelemP2(tt,tj); s13(ind) = Ls13(ti,tj);
                i14(ind) = linelemP2(tt,ti); j14(ind) = elemP2online(tt,tj); s14(ind) = Ls14(ti,tj);
                i15(ind) = linelemP2(tt,ti); j15(ind) = elemP2online(tt,tj); s15(ind) = Ls15(ti,tj);   
                i16(ind) = elemP2online(tt,tj); j16(ind) = linelemP2(tt,ti); s16(ind) = Ls16(ti,tj);
                i17(ind) = elemP2online(tt,tj); j17(ind) = linelemP2(tt,ti); s17(ind) = Ls17(ti,tj);                
            end
        end
    end
    
%% ---------------------------------------------------------------------------
    linelemP2 = [1:Nmem-1; Nmem+1:2*Nmem-1; 2:Nmem]';
    i2 = zeros(3*(Nl-1),1); j2 = zeros(3*(Nl-1),1); s2 = zeros(3*(Nl-1),1);
    ind = 0;
    for t = 1:cl2mem1-1
        len = sqrt((pmem(t+1,1)-pmem(t,1))^2+(pmem(t+1,2)-pmem(t,2))^2);
        Ls2 = len*[kmem(t) kmem(Nmem+t) kmem(t+1)].*[1/6 2/3 1/6];
        for ti = 1:3
            ind = ind + 1;
            i2(ind) = linelemP2(t,ti); j2(ind) = Nm+1; s2(ind) = Ls2(ti);
        end
    end
    

    t = cl2mem1-1;  len = sqrt((pmem(t+1,1)-pmem(t,1))^2+(pmem(t+1,2)-pmem(t,2))^2);
    Ls31 = len/Ca*[kmem(t) kmem(Nmem+t) kmem(t+1)].*[1/6 2/3 1/6].*[nnode_mem(t,1) nelem_mem(t,1) nnode_mem(t+1,1)];
    Ls32 = len/Ca*[kmem(t) kmem(Nmem+t) kmem(t+1)].*[1/6 2/3 1/6].*[nnode_mem(t,2) nelem_mem(t,2) nnode_mem(t+1,2)];
    mem2node = geom.mem2node;
    i31 = mem2node(1); j31 = Nm+1; s31 = Ls31(3);
    i32 = mem2node(1); j32 = Nm+1; s32 = Ls32(3);
    
    stiff_nu = [];
    
    stiff_nu.G11 = sparse(i11,j11,s11,N+Ne,2*Nm);
    stiff_nu.G12 = sparse(i12,j12,s12,N+Ne,2*Nm);
    stiff_nu.G13 = sparse(i13,j13,s13,2*Nm,2*Nm);
    stiff_nu.G14 = sparse(i14,j14,s14,2*Nm,N+Ne);
    stiff_nu.G15 = sparse(i15,j15,s15,2*Nm,N+Ne);
    stiff_nu.G16 = sparse(i16,j16,s16,N+Ne,2*Nm);
    stiff_nu.G17 = sparse(i17,j17,s17,N+Ne,2*Nm);
    
    G2 = sparse(i2,j2,s2,2*Nmem-1,2*Nm);
    stiff_nu.G2 = [G2(1:cl2mem1-1,:); sparse(1,2*Nm); G2(cl2mem1+1:2*Nmem-1,:)];
    
    stiff_nu.G31 = sparse(i31,j31,s31,N+Ne,2*Nm);
    stiff_nu.G32 = sparse(i32,j32,s32,N+Ne,2*Nm);    
    
