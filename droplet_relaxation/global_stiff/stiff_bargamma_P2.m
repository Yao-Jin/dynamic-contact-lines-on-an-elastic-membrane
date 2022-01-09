function [stiff_bargamma] = stiff_bargamma_P2(geom,para,sln)

% bargamma: P2;
% construct stiffness matrices for membrane elements
%% ----------------------------------------------------------------------------
    % read geometry data    
    N = size(geom.node,1); Ne = size(geom.edge,1); Nmem = para.Nmem;
    pmem = sln.pmem(1:Nmem,1:2); mem2node = geom.mem2node; mem2edge = geom.mem2edge;
    cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2);
    
    % precondition constants on each elem
    [~, ~, tauelem_mem, ~] = getnormal_tau_mem_P1(pmem);
    pypxelem = (pmem(2:Nmem,2)-pmem(1:Nmem-1,2))./(pmem(2:Nmem,1)-pmem(1:Nmem-1,1));
    lenelem = sqrt(1+pypxelem.^2);
    
%% ---------------------------------------------------------------------------
    % read model parameters
    ls = para.ls; mu1 = para.mu1; mu2 = para.mu2;

%% ---------------------------------------------------------------------------
    % initialize indices for element
    elemP2online = [mem2node(1:end-1)' mem2edge(1:end)+N mem2node(2:end)'];
    
    i21 = zeros(9*(Nmem-1),1); j21 = zeros(9*(Nmem-1),1); s21 = zeros(9*(Nmem-1),1);
    i22 = zeros(9*(Nmem-1),1); j22 = zeros(9*(Nmem-1),1); s22 = zeros(9*(Nmem-1),1);
    i23 = zeros(9*(Nmem-1),1); j23 = zeros(9*(Nmem-1),1); s23 = zeros(9*(Nmem-1),1);
    i24 = zeros(9*(Nmem-1),1); j24 = zeros(9*(Nmem-1),1); s24 = zeros(9*(Nmem-1),1);  
    
    i6 = zeros(9*(Nmem-1),1); j6 = zeros(9*(Nmem-1),1); s6 = zeros(9*(Nmem-1),1);
    
    linelemP2 = [1:cl2mem1-1; Nmem+3:Nmem+cl2mem1+1; 2:cl2mem1]';
    ind = 0;
    for t = 1:cl2mem1-1
        mu = mu2;
        Ls21 = tauelem_mem(t,1)*getStiffonRef1DP2P2der;
        Ls22 = tauelem_mem(t,2)*getStiffonRef1DP2P2der;
        Ls23 = 1/ls*tauelem_mem(t,1)*getStiffonRef1DP2P2der;
        Ls24 = 1/ls*tauelem_mem(t,2)*getStiffonRef1DP2P2der;
        Ls6 = 1/mu/((pmem(t+1,1)-pmem(t,1))*lenelem(t))*getStiffonRef1DP2der;
        for ti = 1:3
            for tj = 1:3
                ind = ind + 1;
                i21(ind) = linelemP2(t,ti); j21(ind) = elemP2online(t,tj); s21(ind) = Ls21(ti,tj);
                i22(ind) = linelemP2(t,ti); j22(ind) = elemP2online(t,tj); s22(ind) = Ls22(ti,tj);
                i23(ind) = linelemP2(t,ti); j23(ind) = elemP2online(t,tj); s23(ind) = Ls23(ti,tj);
                i24(ind) = linelemP2(t,ti); j24(ind) = elemP2online(t,tj); s24(ind) = Ls24(ti,tj); 
                i6(ind) = linelemP2(t,ti); j6(ind) = linelemP2(t,tj); s6(ind) = Ls6(ti,tj);  
            end
        end
    end
    
    linelemP2 = [cl2mem2:Nmem-1; Nmem+cl2mem2+2:2*Nmem+1; cl2mem2+1:Nmem]';
    for t = cl2mem2:Nmem-1
        mu = mu2;
        Ls21 = tauelem_mem(t,1)*getStiffonRef1DP2P2der;
        Ls22 = tauelem_mem(t,2)*getStiffonRef1DP2P2der;
        Ls23 = 1/ls*tauelem_mem(t,1)*getStiffonRef1DP2P2der;
        Ls24 = 1/ls*tauelem_mem(t,2)*getStiffonRef1DP2P2der;       
        Ls6 = 1/mu/((pmem(t+1,1)-pmem(t,1))*lenelem(t))*getStiffonRef1DP2der;
        for ti = 1:3
            for tj = 1:3
                ind = ind + 1;
                i21(ind) = linelemP2(t-cl2mem2+1,ti); j21(ind) = elemP2online(t,tj); s21(ind) = Ls21(ti,tj);
                i22(ind) = linelemP2(t-cl2mem2+1,ti); j22(ind) = elemP2online(t,tj); s22(ind) = Ls22(ti,tj);
                i23(ind) = linelemP2(t-cl2mem2+1,ti); j23(ind) = elemP2online(t,tj); s23(ind) = Ls23(ti,tj);
                i24(ind) = linelemP2(t-cl2mem2+1,ti); j24(ind) = elemP2online(t,tj); s24(ind) = Ls24(ti,tj);                               
                i6(ind) = linelemP2(t-cl2mem2+1,ti); j6(ind) = linelemP2(t-cl2mem2+1,tj); s6(ind) = Ls6(ti,tj);  
            end
        end
    end
    
    linelemP2 = [Nmem+1 cl2mem1+1:cl2mem2-1; Nmem+cl2mem1+2:Nmem+cl2mem2+1; cl2mem1+1:cl2mem2-1 Nmem+2]';
    for t = cl2mem1:cl2mem2-1
        mu = mu1;
        Ls21 = tauelem_mem(t,1)*getStiffonRef1DP2P2der;
        Ls22 = tauelem_mem(t,2)*getStiffonRef1DP2P2der;
        Ls23 = 1/ls*tauelem_mem(t,1)*getStiffonRef1DP2P2der;
        Ls24 = 1/ls*tauelem_mem(t,2)*getStiffonRef1DP2P2der;     
        Ls6 = 1/mu/((pmem(t+1,1)-pmem(t,1))*lenelem(t))*getStiffonRef1DP2der;
        for ti = 1:3
            for tj = 1:3
                ind = ind + 1;
                i21(ind) = linelemP2(t-cl2mem1+1,ti); j21(ind) = elemP2online(t,tj); s21(ind) = Ls21(ti,tj);
                i22(ind) = linelemP2(t-cl2mem1+1,ti); j22(ind) = elemP2online(t,tj); s22(ind) = Ls22(ti,tj);
                i23(ind) = linelemP2(t-cl2mem1+1,ti); j23(ind) = elemP2online(t,tj); s23(ind) = Ls23(ti,tj);
                i24(ind) = linelemP2(t-cl2mem1+1,ti); j24(ind) = elemP2online(t,tj); s24(ind) = Ls24(ti,tj);                  
                i6(ind) = linelemP2(t-cl2mem1+1,ti); j6(ind) = linelemP2(t-cl2mem1+1,tj); s6(ind) = Ls6(ti,tj);  
            end
        end
    end
    
    stiff_bargamma = [];
    
    stiff_bargamma.M21 = sparse(i21,j21,s21,2*Nmem+1,N+Ne);
    stiff_bargamma.M22 = sparse(i22,j22,s22,2*Nmem+1,N+Ne);
    
    stiff_bargamma.M23 = sparse(i23,j23,s23,2*Nmem+1,N+Ne);
    stiff_bargamma.M24 = sparse(i24,j24,s24,2*Nmem+1,N+Ne);
    
    stiff_bargamma.M6 = sparse(i6,j6,s6,2*Nmem+1,2*Nmem+1);
    
 