function [geom] = Movingmesh_lrperiodic(geom,para,sln)
% moving mesh method
% construct the new mesh based on the old one but keep the mesh connectivity and topology

    elem = geom.elem; node = geom.node;
    N = size(node,1); Nt = size(elem,1);
    Nmem = para.Nmem; Nffi = para.Nffi;
    
% calculate the specific parameter for the moving mesh method
    p1 = node(elem(1:Nt,1),1:2);
    p2 = node(elem(1:Nt,2),1:2);
    p3 = node(elem(1:Nt,3),1:2);
    area = 0.5*abs((p2(:,1)-p1(:,1)).*(p3(:,2)-p1(:,2))-(p3(:,1)-p1(:,1)).*(p2(:,2)-p1(:,2)));
    maxarea = max(area); minarea = min(area);
    sp = 1+(maxarea-minarea)./area;
    
    pmem = sln.pmem(1:Nmem,1:2); pffi = sln.pffi(1:Nffi,1:2);
    ffi2node = geom.ffi2node; mem2node = geom.mem2node;
    lbd2node = geom.lbd2node; rbd2node = geom.rbd2node;
    bdnodeup = geom.bdnodeup; 
    
% initialize the solution including boundary conditions
    solved = zeros(N,2);
    solved(mem2node,1:2) = pmem - sln.pmemold(1:Nmem,1:2);
    solved(ffi2node,1:2) = pffi - sln.pffiold(1:Nffi,1:2);
    
% construct stiffness matrix based on each element

    i11 = zeros(9*Nt,1); j11 = zeros(9*Nt,1); s11 = zeros(9*Nt,1);
    i12 = zeros(9*Nt,1); j12 = zeros(9*Nt,1); s12 = zeros(9*Nt,1);
    i13 = zeros(9*Nt,1); j13 = zeros(9*Nt,1); s13 = zeros(9*Nt,1);
    i14 = zeros(9*Nt,1); j14 = zeros(9*Nt,1); s14 = zeros(9*Nt,1);
    i21 = zeros(9*Nt,1); j21 = zeros(9*Nt,1); s21 = zeros(9*Nt,1);
    i22 = zeros(9*Nt,1); j22 = zeros(9*Nt,1); s22 = zeros(9*Nt,1);
    i23 = zeros(9*Nt,1); j23 = zeros(9*Nt,1); s23 = zeros(9*Nt,1);
    i24 = zeros(9*Nt,1); j24 = zeros(9*Nt,1); s24 = zeros(9*Nt,1);

    ind = 0;
    for t = 1:Nt
        [Ls11, Ls12, Ls13, Ls14] = localstiff_movmesh1(node(elem(t,:),:),sp(t));
        [Ls21, Ls22, Ls23, Ls24] = localstiff_movmesh2(node(elem(t,:),:),sp(t));
        for i = 1:3
            for j = 1:3
                ind = ind+1;
                i11(ind) = elem(t,i); j11(ind) = elem(t,j); s11(ind) = Ls11(i,j);
                i12(ind) = elem(t,i); j12(ind) = elem(t,j); s12(ind) = Ls12(i,j);
                i13(ind) = elem(t,i); j13(ind) = elem(t,j); s13(ind) = Ls13(i,j);
                i14(ind) = elem(t,i); j14(ind) = elem(t,j); s14(ind) = Ls14(i,j);
                i21(ind) = elem(t,i); j21(ind) = elem(t,j); s21(ind) = Ls21(i,j);
                i22(ind) = elem(t,i); j22(ind) = elem(t,j); s22(ind) = Ls22(i,j);
                i23(ind) = elem(t,i); j23(ind) = elem(t,j); s23(ind) = Ls23(i,j);
                i24(ind) = elem(t,i); j24(ind) = elem(t,j); s24(ind) = Ls24(i,j);
            end
        end
    end
    
    A11 = sparse(i11,j11,s11,N,N);
    A12 = sparse(i12,j12,s12,N,N);
    A13 = sparse(i13,j13,s13,N,N);
    A14 = sparse(i14,j14,s14,N,N);
    A21 = sparse(i21,j21,s21,N,N);
    A22 = sparse(i22,j22,s22,N,N);
    A23 = sparse(i23,j23,s23,N,N);
    A24 = sparse(i24,j24,s24,N,N);

    
    tempmat = [A11+A21 A12+A22;
               A13+A23 A14+A24];
    
    modifynode = ones(N,1);
    modifynode(ffi2node) = 0; modifynode(mem2node) = 0;
    modifynode(bdnodeup(:)==1) = 0;
    modifynodey = modifynode;
    modifynode(lbd2node) = 0; modifynode(rbd2node) = 0;
    modifynodex = modifynode;
    
    indexkeep = [modifynodex;modifynodey];       
    
    tempload = -tempmat(:,indexkeep(:)==0)*[solved(modifynodex(:)==0,1);solved(modifynodey(:)==0,2)];       
           
    indexl = [lbd2node;lbd2node+N];
    indexr = [rbd2node;rbd2node+N];
    
    tempmat(indexl,:) = tempmat(indexl,:)+tempmat(indexr,:);
    tempmat(:,indexl) = tempmat(:,indexl)+tempmat(:,indexr);
    tempload(indexl) = tempload(indexl)+tempload(indexr);
    
    modifynodey(rbd2node) = 0;
    
    indexkeep = [modifynodex;modifynodey];
    stiffmat = tempmat(indexkeep(:)==1,indexkeep(:)==1);
    load = tempload(indexkeep(:)==1);
       
    tempsolved = stiffmat\load;
    solved(modifynodex(:)==1,1) = tempsolved(1:sum(modifynodex));
    solved(modifynodey(:)==1,2) = tempsolved(sum(modifynodex)+1:end);
    
    solved(rbd2node,2) = solved(lbd2node,2);
    
    newnode = node+solved; 
    geom.node = newnode;

    
%     trimesh(geom.elem,newnode(:,1),newnode(:,2),zeros(N,1));
%     hold on;
%     plot(node(:,1),node(:,2),'or');
%     view(2);
%     axis equal;
%     hold off;
%     pause(0.001);
    
end
    