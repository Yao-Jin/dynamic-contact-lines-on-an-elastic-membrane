function [geom] = produce_newmesh(geom,para,sln)
    
    addpath(genpath('.'));
%---------------------------------------------------------------------------------
    Nffi = para.Nffi; Nmem = para.Nmem;
    pffi = sln.pffi(1:Nffi,:); pmem = sln.pmem(1:Nmem,:);
    cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2);
    Nlrbd = para.Nlrbd; 
   
    lbd2node = geom.lbd2node; rbd2node = geom.rbd2node;
    leftbd = geom.node(lbd2node,:); rightbd = geom.node(rbd2node,:);
    hmem = 0;
    hmem = max(hmem, rightbd(2,2)-rightbd(1,2));
    len = sqrt((pffi(2:end,1)-pffi(1:end-1,1)).^2+((pffi(2:end,2)-pffi(1:end-1,2)).^2));
    hmem = max(hmem, max(len));  
    len = sqrt((pmem(2:end,1)-pmem(1:end-1,1)).^2+((pmem(2:end,2)-pmem(1:end-1,2)).^2));
    hmem = max(hmem, max(len));      
    
    polybd = [pmem(1:Nmem,:);1 1;-1 1];
    pfix = [pffi; pmem(1:cl2mem1-1,:); pmem(cl2mem1+1:cl2mem2-1,:); pmem(cl2mem2+1:end,:); -1 1; 1 1; leftbd; rightbd];
    adapt_line = [pmem(1:cl2mem1-1,:);pffi;pmem(cl2mem2+1:Nmem,:);
                  pmem(Nmem:-1:1,:)+[zeros(Nmem,1) -0.0001*ones(Nmem,1)]];
              
    disp('reproducing mesh......');
    Bx = para.Bx; By = para.By;
    find_in_out = true; fd = @(p)(p_poly_dist(p,polybd,find_in_out));
    fhout = @(p)(min(hmem*(1+20*abs(p_poly_dist(p,adapt_line)).^2),5*hmem));
%     fhout = @(p)(min(hmem*(1+16*abs(p_poly_dist(p,adapt_line)).^2),4*hmem));
%     fhout = @(p)(hmem*(1+4*abs(p_poly_dist(p,adapt_line))));
%     fhout = @(p)(hmem*(1+6*abs(p_poly_dist(p,adapt_line)).^2));
    [node,elem] = mydistmesh2d(fd,fhout,hmem,[-Bx, min(polybd(:,2));Bx,By],pfix);
    
    [geom] = auxstructure(node,elem,para.Bx,para.By);
    geom.cl2mem = [cl2mem1 cl2mem2];
    
    %----------------------------------------------------------------------------
    trimesh(elem,node(:,1),node(:,2),zeros(size(node,1),1));
    hold on;
    plot(pfix(:,1),pfix(:,2),'r.');
    plot(node(geom.bdnode(:)==1,1),node(geom.bdnode(:)==1,2),'bo');
    view(2);
    axis equal;
    hold off;
    %----------------------------------------------------------------------------
    
    % identify nodes on ffi
    ffi2node = (1:Nffi);
    geom.ffi2node = ffi2node;
    geom.cl2node = [1 Nffi];
    
    % identify nodes on membrane
    mem2node = [Nffi+1:Nffi+cl2mem1-1 1 Nffi+cl2mem1:Nffi+cl2mem2-2 Nffi Nffi+cl2mem2-1:Nffi+Nmem-2];
    geom.mem2node = mem2node;

    % identify nodes on the left boundary
    lbd2node= (Nffi+Nmem+1:Nffi+Nmem+Nlrbd)';
    geom.lbd2node = lbd2node;

    % identify nodes on the right boundary
    rbd2node= (Nffi+Nmem+Nlrbd+1:Nffi+Nmem+2*Nlrbd)';
    geom.rbd2node = rbd2node;

    % identify nodes on the upper boundary
    upbd2node = find(geom.bdnodeup(:)==1);
    para.Nupbd = length(upbd2node);
    tmp = geom.node(upbd2node,1);
    for i=1:para.Nupbd
        for j=i+1:para.Nupbd
            if tmp(i)>tmp(j)
                stand=tmp(i); tmp(i)=tmp(j); tmp(j)=stand;
                stand=upbd2node(i); upbd2node(i)=upbd2node(j); upbd2node(j)=stand;
            end
        end
    end
    geom.upbd2node = upbd2node;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    edge = geom.edge;
    e2v = geom.e2v;
    % identify edges for ffi
    edge2ffi = zeros(size(edge,1),1);
    ffi2edge = zeros(Nffi-1,1);
    for i = 1:Nffi-1
        edge2ffi(e2v(:,ffi2node(i)).*e2v(:,ffi2node(i+1))==1)=1;
        ffi2edge(i) = find(e2v(:,ffi2node(i)).*e2v(:,ffi2node(i+1))==1);
    end

    geom.edgeffi = edge2ffi;
    geom.ffi2edge = ffi2edge;

    % identify edges for mem
    edge2mem = geom.bdedgedown;
    geom.edgemem = edge2mem;
    mem2edge = zeros(Nmem-1,1);
    for i = 1:Nmem-1
        mem2edge(i) = find(e2v(:,mem2node(i)).*e2v(:,mem2node(i+1))==1);
    end
    geom.mem2edge = mem2edge;

    % identify edges on the left boundary
    lbd2node=[mem2node(1); lbd2node; upbd2node(1)];
    lbd2edge = zeros(Nlrbd+1,1);
    for i = 1:Nlrbd+1
        lbd2edge(i) = find(e2v(:,lbd2node(i)).*e2v(:,lbd2node(i+1))==1);
    end
    geom.lbd2edge = lbd2edge;

    % identify edges on the right boundary
    rbd2node=[mem2node(end); rbd2node; upbd2node(end)];
    rbd2edge = zeros(Nlrbd+1,1);
    for i = 1:Nlrbd+1
        rbd2edge(i) = find(e2v(:,rbd2node(i)).*e2v(:,rbd2node(i+1))==1);
    end
    geom.rbd2edge = rbd2edge;
    
    % identify edges on the upper boundary
    upbd2edge = zeros(para.Nupbd-1,1);
    for i = 1:para.Nupbd-1
        upbd2edge(i) = find(e2v(:,upbd2node(i)).*e2v(:,upbd2node(i+1))==1);
    end
    geom.upbd2edge = upbd2edge;

    % judge elem inside or outside the droplet
    bd_droplet = [pmem(cl2mem1+1:cl2mem2-1,1:2); pffi(Nffi:-1:1,1:2)];
    [inout] = judge_inout(node,elem,bd_droplet);
    geom.inout = inout;
    
end