function [stiff_bulk] = stiff_Bulk_P2P01(geom,para,sln)

% construct stiffness matrices for bulk elements

%% ----------------------------------------------------------------------------
    % read geometry data
    elem = geom.elem; node = geom.node; edge = geom.edge;
    Nt = size(elem,1); N = size(node,1); Ne = size(edge,1);
    elem2edge = geom.elem2edge;
    
    cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2);
    node(geom.mem2node,2) = sln.pmem(cl2mem1:cl2mem2,2);
%% ---------------------------------------------------------------------------
    % read model parameters
    eta = para.eta1;
    
%% ---------------------------------------------------------------------------
    % initialize indices for element 
    elemP2 = [elem elem2edge(:,3)+N elem2edge(:,1)+N elem2edge(:,2)+N]; 
%     elemP0 = linspace(1,Nt,Nt)';
    elemP01 = [linspace(1,Nt,Nt)'+N elem];

    %% based on triangular elements 
    i11 = zeros(36*Nt,1); j11 = zeros(36*Nt,1); s11 = zeros(36*Nt,1);
    i12 = zeros(36*Nt,1); j12 = zeros(36*Nt,1); s12 = zeros(36*Nt,1);
    i13 = zeros(36*Nt,1); j13 = zeros(36*Nt,1); s13 = zeros(36*Nt,1);
    i14 = zeros(36*Nt,1); j14 = zeros(36*Nt,1); s14 = zeros(36*Nt,1);
    i21 = zeros(24*Nt,1); j21 = zeros(24*Nt,1); s21 = zeros(24*Nt,1);
    i22 = zeros(24*Nt,1); j22 = zeros(24*Nt,1); s22 = zeros(24*Nt,1);
    
    ind1 = 0; ind2 = 0;
    for t = 1:Nt
        
        [Ls11, Ls12, Ls13, Ls14] = localstiff_bulk1(node(elem(t,:),:),eta);
        [Ls21, Ls22] = localstiff_bulk2(node(elem(t,:),:));
        
        for ti = 1:6
            for tj = 1:6
                ind1 = ind1 + 1;
                i11(ind1) = elemP2(t,ti); j11(ind1) = elemP2(t,tj); s11(ind1) = Ls11(ti,tj);
                i12(ind1) = elemP2(t,ti); j12(ind1) = elemP2(t,tj); s12(ind1) = Ls12(ti,tj);
                i13(ind1) = elemP2(t,ti); j13(ind1) = elemP2(t,tj); s13(ind1) = Ls13(ti,tj);
                i14(ind1) = elemP2(t,ti); j14(ind1) = elemP2(t,tj); s14(ind1) = Ls14(ti,tj);
            end
        end
    
        for ti = 1:4
            for tj = 1:6
                ind2 = ind2 + 1;
                i21(ind2) = elemP01(t,ti); j21(ind2) = elemP2(t,tj); s21(ind2) = Ls21(ti,tj);
                i22(ind2) = elemP01(t,ti); j22(ind2) = elemP2(t,tj); s22(ind2) = Ls22(ti,tj);
            end
        end
        
    end
    
    stiff_bulk = [];
    stiff_bulk.B11 = sparse(i11,j11,s11,N+Ne,N+Ne);
    stiff_bulk.B12 = sparse(i12,j12,s12,N+Ne,N+Ne);
    stiff_bulk.B13 = sparse(i13,j13,s13,N+Ne,N+Ne);
    stiff_bulk.B14 = sparse(i14,j14,s14,N+Ne,N+Ne);
    stiff_bulk.B21 = sparse(i21,j21,s21,Nt+N,N+Ne);
    stiff_bulk.B22 = sparse(i22,j22,s22,Nt+N,N+Ne);
        
end
