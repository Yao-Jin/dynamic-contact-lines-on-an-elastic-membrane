function [sln, para] = fluid_P2P01_mem_P2P2_bargamma_P2(geom,para,sln)

%  construct and solve the coupled system

    elem = geom.elem; node = geom.node; edge = geom.edge; 
    Nt = size(elem,1); N = size(node,1); Ne = size(edge,1); Nmem = para.Nmem; Nffi = para.Nffi;
    cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2); 
    pmem = sln.pmem(:,1:2); pffi = sln.pffi(:,1:2); dt = para.dt;
    
% construct the stiffness matrix and the load vector by assembling local stiffness matrices

% stiffness matrix for bulk elements
    [stiff_bulk] = stiff_Bulk_P2P01(geom,para,sln);
    B11 = stiff_bulk.B11; B12 = stiff_bulk.B12; B13 = stiff_bulk.B13; B14 = stiff_bulk.B14;
    B21 = stiff_bulk.B21; B22 = stiff_bulk.B22;
    
% stiffness matrix for fluid interface elements
    [stiff_interface] = stiff_Interface_P2(geom,para,sln);
    I1 = stiff_interface.I1;

% stiffness matrix for membrane elements
    [stiff_membrane] = stiff_Membrane_P2P2(geom,para,sln);
    M11 = stiff_membrane.M11; M12 = stiff_membrane.M12;
    M31 = stiff_membrane.M31; M32 = stiff_membrane.M32; 
    M41 = stiff_membrane.M41; M42 = stiff_membrane.M42; 
    M43 = stiff_membrane.M43; M44 = stiff_membrane.M44;
    M45 = stiff_membrane.M45;
    M51 = stiff_membrane.M51; M52 = stiff_membrane.M52;
 
    [stiff_bargamma] = stiff_bargamma_P2(geom,para,sln);
    M21 = stiff_bargamma.M21; M22 = stiff_bargamma.M22; 
    M23 = stiff_bargamma.M23; M24 = stiff_bargamma.M24; M6 = stiff_bargamma.M6;    
    
% stiffness matrix for elements at the contact line  
    [stiff_cl] = stiff_Contactline(geom,para,sln);
    cl1 = stiff_cl.cl1; cl4 = stiff_cl.cl4; 
    
    % u_x u_y p y kappa bargamma_mu r_x r_y
    % here we are solving the difference of y and r_x, r_y from time t^m to t^(m+1)
    N1 = 2*Nmem-1; N3 = Nt+N; N4 = N+Ne; N5 = 2*Nmem+1;
    
    tempStiff = [B11 B12 -B21' -M51 M31 -M21';
                 B13 B14 -B22' -M52 M32 -M22';
                 B21 B22 sparse(N3,N3+2*N1+N5);
                 M42 M43 sparse(N1,N3) M41/dt sparse(N1,N1+N5);
                 sparse(N1,2*N4+N3) M45 M44 sparse(N1,N5);
                 M23 M24 sparse(N5,N3+2*N1) M6+cl4];
             
    tempload = [-cl1;sparse(N4+N3+2*N1+N5,1)];
    
    tempload(1:N4) = tempload(1:N4)-I1*pffi(:,1)+(-M11+M51)*pmem(:,2);
    tempload(N4+1:2*N4) = tempload(N4+1:2*N4)-I1*pffi(:,2)+(-M12+M52)*pmem(:,2);
    tempload(2*N4+N3+N1+1:2*N4+N3+2*N1) = tempload(2*N4+N3+N1+1:2*N4+N3+2*N1)-M45*pmem(:,2);    
    
%% --------------------------------------------------------------------------------
    % extract unknown variables' indices
 
    % bdc of fluid
    bdnodeup = geom.bdnodeup; bdedgeup = geom.bdedgeup;
    
    ind1x = [bdnodeup;bdedgeup]; ind1x(geom.mem2node(floor((Nmem+1)/2))) = 1;
    ind1y = [bdnodeup;bdedgeup];
    ind1p = zeros(N3,1); ind1p(1) = 1;
    
    ind2F = [0;zeros(Nmem-2,1);0;zeros(Nmem-1,1)]; ind2D = [1;zeros(Nmem-2,1);1;zeros(Nmem-1,1)]; 
    ind2D2 = [1;zeros(Nmem-2,1);1;0;0;zeros(Nmem-1,1)];

    indexdelete = [ind1x;ind1y;ind1p;ind2F;ind2D;ind2D2];

    Stiffmat = tempStiff(indexdelete(:)==0, indexdelete(:)==0);

    load = tempload(indexdelete(:)==0);
    
    solved = Stiffmat\load;
    
%% -----------------------------------------------------------------------------
    % arrange solutions
    back = zeros(size(indexdelete));
    back(indexdelete(:)==0)=solved;

    sln.pmemold = sln.pmem; sln.kmemold = sln.kmem; 
    sln.bargamma_old = sln.bargamma; sln.pffiold = sln.pffi;
    
    index = 2*N4+N3;
    sln.pmem(:,2) = pmem(:,2)+back(index+1:index+N1); 
    sln.maxdy = max(abs(back(index+1:index+N1)));
    sln.dydt = back(index+1:index+N1)/dt;
    sln.kmem = back(index+N1+1:index+2*N1);
    sln.bargamma = back(index+2*N1+1:index+2*N1+N5);    
    
    % update the contact lines' height
    sln.pffi(1,2) = sln.pmem(cl2mem1,2); sln.pffi(Nffi,2) = sln.pmem(cl2mem2,2);
    
    flu_vel = [back(1:N4) back(N4+1:2*N4)];
    sln.flu_vel = flu_vel;
    pressure = back(2*N4+1:2*N4+N3);
    sln.pressure = pressure;
   
    % visualization
    visual = false;
    while visual
        pos = [node; 0.5*(node(geom.edge(:,1),1:2)+node(geom.edge(:,2),1:2))];
        trimesh(elem,node(:,1),node(:,2),zeros(size(node,1),1));
        hold on;
        quiver(pos(:,1),pos(:,2),flu_vel(:,1),flu_vel(:,2),'r');
        plot(pmem(cl2mem1,1),pmem(cl2mem1,2),'bo');
        view(2);
        axis equal;
        hold off;
        pause(0.01);
               
        elem_mid(:,1) = 1/3*(node(elem(:,1),1)+node(elem(:,2),1)+node(elem(:,3),1));
        elem_mid(:,2) = 1/3*(node(elem(:,1),2)+node(elem(:,2),2)+node(elem(:,3),2));
        plot3(elem_mid(:,1),elem_mid(:,2),pressure(N+1:N+Nt),'or'); 
        
        pmem_nodal = sln.pmem(1:Nmem,:);
        pffi_nodal = sln.pffi(1:Nffi,:);
        
        pffi_mid(:,1) = 1/2*(pffi_nodal(1:end-1,1)+pffi_nodal(2:end,1));
        pffi_mid(:,2) = 1/2*(pffi_nodal(1:end-1,2)+pffi_nodal(2:end,2));
        pmem_mid(:,1) = 1/2*(pmem_nodal(1:end-1,1)+pmem_nodal(2:end,1));
        pmem_mid(:,2) = 1/2*(pmem_nodal(1:end-1,2)+pmem_nodal(2:end,2));
 
        %--------------------------------------------------------------
        %check the symmetry
%         disp(['pmem symmetry check px: ',num2str(max(abs(pmem_nodal(1:Nmem,1)+pmem_nodal(Nmem:-1:1,1))))]);
%         disp(['pmem symmetry check py: ',num2str(max(abs(pmem_nodal(1:Nmem,2)-pmem_nodal(Nmem:-1:1,2))))]);
%         disp(['kmem symmetry check: ',num2str(max(abs(sln.kmem(1:Nmem)-sln.kmem(Nmem:-1:1))))]);
%         disp(['bargamma symmetry check: ',num2str(max(abs(sln.bargamma(1:Nmem)-sln.bargamma(Nmem:-1:1))))]);
        %---------------------------------------------------------------
        
        plot(pmem_nodal(:,1),pmem_nodal(:,2),'o-r');
        hold on;
        plot(pffi_nodal(:,1),pffi_nodal(:,2),'o-b');
        plot(node(1,1),node(1,2),'g*');
        plot(node(Nffi,1),node(Nffi,2),'g*');
        axis tight;
        ffi2edge = geom.ffi2edge; mem2edge = geom.mem2edge;
        quiver([node(:,1);pffi_mid(:,1);pmem_mid(:,1)],[node(:,2);pffi_mid(:,2);pmem_mid(:,2)], ...
           [flu_vel(1:N,1);flu_vel(ffi2edge+N,1);flu_vel(mem2edge+N,1)],...
           [flu_vel(1:N,2);flu_vel(ffi2edge+N,2);flu_vel(mem2edge+N,2)],'k');
        hold off;
%         pause(0.01);
        
        plot(pmem_nodal(:,1),pmem_nodal(:,2),'o-r');
        hold on;
        plot(pffi_nodal(:,1),pffi_nodal(:,2),'o-r');
        plot(sln.pmemold(1:Nmem,1),sln.pmemold(1:Nmem,2),'*--b');
        plot(sln.pffiold(1:Nffi,1),sln.pffiold(1:Nffi,2),'*--b');
        axis equal;
        hold off;
        
        plot(pmem_nodal(:,1),sln.kmem(1:Nmem),'-*r');
        hold on;
        plot(sln.pmem(Nmem+1:end,1),sln.kmem(Nmem+1:end),'-sb');
        hold off;
%         pause(0.01);
        
        pmem_nodal = sln.pmemold(1:Nmem,1:2);
        plot(pmem_nodal(1:cl2mem1,1),sln.bargamma(1:cl2mem1),'-ob');
        hold on;
        plot(pmem_nodal(cl2mem1:cl2mem2,1),[sln.bargamma(Nmem+1);sln.bargamma(cl2mem1+1:cl2mem2-1);sln.bargamma(Nmem+2)],'-sr');
        plot(pmem_nodal(cl2mem2:Nmem,1),sln.bargamma(cl2mem2:Nmem),'-ob');    
        hold off;     
%         pause(0.001);
        
%         disp(['left bargamma jump: ', num2str(sln.bargamma(Nmem+1)-sln.bargamma(cl2mem1))]);    
%         disp(['right bargamma jump: ', num2str(sln.bargamma(Nmem+2)-sln.bargamma(cl2mem2))]); 

        break;
    end
     
     para.curlen = getlength(sln.pmem(1:Nmem,:));
%      disp(['curlen: ',num2str(para.curlen)]); 
     bd_droplet = [sln.pffi(1:Nffi,:); sln.pmem(cl2mem2-1:-1:cl2mem1+1,:)];
     areanow = polyarea(bd_droplet(:,1),bd_droplet(:,2));
     para.areanow = areanow;
     disp(['droplet area after mem: ',num2str(areanow)]);     

end