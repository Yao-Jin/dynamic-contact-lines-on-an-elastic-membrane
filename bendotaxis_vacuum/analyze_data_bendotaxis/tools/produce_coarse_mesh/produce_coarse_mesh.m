function [node,elem] = produce_causemesh(geom,para,sln)
    
%---------------------------------------------------------------------------------
    Nffil = para.Nffil; Nffir = para.Nffir;
    pffil = sln.pffilold; pffir = sln.pffirold; pmem = sln.pmemold;
    cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2); 
    By = para.By;
    
    Nc = 5;
    h = getlength(pffir(1:Nffir,1:2))/Nc;
    
    cpffir = interparc(Nc,pffir(1:Nffir,1),pffir(1:Nffir,2));
    cpffil = interparc(round(getlength(pffil(1:Nffil,1:2))/h),pffil(1:Nffil,1),pffil(1:Nffil,2));
    cpmem = interparc(round(getlength(pmem(cl2mem1:cl2mem2,1:2))/h),pmem(cl2mem1:cl2mem2,1),pmem(cl2mem1:cl2mem2,2));
    
    polybd = [cpffil; cpffir(end:-1:1,:); cpmem(end-1:-1:2,1:2)];
    pfix = [cpffil; cpffir; cpmem(2:end-1,1:2)];             
%     adapt_line = [cpffil(:,1)-0.00001 cpffil(:,2); cpffil(end:-1:1,1)+0.00001 cpffil(end:-1:1,2);
%                   cpmem(2:end-1,:);
%                   cpffir(:,1)-0.00001 cpffir(:,2); cpffir(end:-1:1,1)+0.00001 cpffir(end:-1:1,2);
%                   cpmem(end-1:-1:2,1) cpmem(end-1:-1:2,2)-0.00001];
    adapt_line = [cpmem; cpmem(end:-1:1,1) cpmem(end:-1:1,2)-0.00001];
    disp('reproducing mesh......');
    fd = @(p)(p_poly_dist(p,polybd));
    fhout = @(p)(h*(1+0*abs(p_poly_dist(p,adapt_line))));
    [node,elem] = mydistmesh2d(fd,fhout,h,[min(cpffil(:,1)),min(cpmem(:,2));max(cpffir(:,1)),By],pfix);
    
    %----------------------------------------------------------------------------
    trimesh(elem,node(:,1),node(:,2),zeros(size(node,1),1));
    hold on;
    plot(pfix(:,1),pfix(:,2),'r.');
    view(2);
    axis equal;
    hold off;
    %----------------------------------------------------------------------------
    
end