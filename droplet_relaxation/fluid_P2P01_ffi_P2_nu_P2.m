function [sln, para] = fluid_P2P01_ffi_P2_nu_P2(geom,para,sln)

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
    I1 = stiff_interface.I1; I2 = stiff_interface.I2; I3 = stiff_interface.I3;
    
% stiffness matrix for membrane elements
    [stiff_nu] = stiff_nu_P2(geom,para,sln);
    M21 = stiff_nu.M21; M22 = stiff_nu.M22;
    M23 = stiff_nu.M23; M24 = stiff_nu.M24; M6 = stiff_nu.M6; 
    
% stiffness matrix for elements at the contact line  
    [stiff_cl] = stiff_Contactline(geom,para,sln);
    cl4 = stiff_cl.cl4; cl5 = stiff_cl.cl5; cl6 = stiff_cl.cl6; 
    
    % u_x u_y p y kappa nu r_x r_y
    % here we are solving the difference of y and r_x, r_y from time t^m to t^(m+1)
    N2 = 2*Nffi-1; N3 = Nt+N; N4 = N+Ne; N5 = 2*Nmem+1;
    
    tempStiff = [B11 B12 -B21' I1 sparse(N4,N2) -M21';
                 B13 B14 -B22' sparse(N4,N2) I1 -M22';
                 B21 B22 sparse(N3,N3+2*N2+N5);
                 -I3' sparse(N2,N4+N3) I2/dt sparse(N2,N2+N5);
                 sparse(N2,N4) -I3' sparse(N2,N3+N2) I2/dt sparse(N2,N5);
                 M23 M24 sparse(N5,N3+2*N2) M6+cl4];    
   
% contact points move along the membrane             
    Dy1 = 1./(pmem(cl2mem1+1,1)-pmem(cl2mem1-1,1))*[pmem(cl2mem1+1,1)-pmem(cl2mem1,1) pmem(cl2mem1,1)-pmem(cl2mem1-1,1)]*...
          [(pmem(cl2mem1,2)-pmem(cl2mem1-1,2))/(pmem(cl2mem1,1)-pmem(cl2mem1-1,1));(pmem(cl2mem1+1,2)-pmem(cl2mem1,2))/(pmem(cl2mem1+1,1)-pmem(cl2mem1,1))];

    Dy2 = 1./(pmem(cl2mem2+1,1)-pmem(cl2mem2-1,1))*[pmem(cl2mem2+1,1)-pmem(cl2mem2,1) pmem(cl2mem2,1)-pmem(cl2mem2-1,1)]*...
          [(pmem(cl2mem2,2)-pmem(cl2mem2-1,2))/(pmem(cl2mem2,1)-pmem(cl2mem2-1,1));(pmem(cl2mem2+1,2)-pmem(cl2mem2,2))/(pmem(cl2mem2+1,1)-pmem(cl2mem2,1))];
       
    tempStiff(2*N4+N3+N2+1,:) = [sparse(1,2*N4+N3) -Dy1/dt sparse(1,N2-1) 1/dt sparse(1,N2-1) sparse(1,N5)];
    tempStiff(2*N4+N3+N2+Nffi,:) = [sparse(1,2*N4+N3) sparse(1,Nffi-1) -Dy2/dt sparse(1,Nffi-1) sparse(1,Nffi-1) 1/dt sparse(1,Nffi-1) sparse(1,N5)];
    
    tempload = [-cl5;-cl6;sparse(N3+2*N2+N5,1)];
       
    tempload(1:N4) = tempload(1:N4)-I1*pffi(:,1);
    tempload(N4+1:2*N4) = tempload(N4+1:2*N4)-I1*pffi(:,2);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Dirichlet boundary condition of fluid velocity along the membrane at the normal direction
    [nelem, nnode, tauelem, taunode] = getnormal_tau_mem_P1(pmem(1:Nmem,1:2));
    dydt = sln.dydt; mem2node = geom.mem2node; mem2edge = geom.mem2edge;
    
    tempStiff(mem2node,:) = taunode(1:Nmem,1).*tempStiff(mem2node,:)...
                          + taunode(1:Nmem,2).*tempStiff(N4+mem2node,:);
    tempload(mem2node) = taunode(1:Nmem,1).*tempload(mem2node)+taunode(1:Nmem,2).*tempload(N4+mem2node);             
    
    tempStiff(N4+mem2node,:) = sparse([1:Nmem 1:Nmem]',[mem2node N4+mem2node]', [nnode(:,1); nnode(:,2)],Nmem,2*N4+N3+2*N2+N5);
    tempload(N4+mem2node) = dydt(1:Nmem).*nnode(:,2);    
    
    tempStiff(N+mem2edge,:) = tauelem(1:Nmem-1,1).*tempStiff(N+mem2edge,:)...
                               + tauelem(1:Nmem-1,2).*tempStiff(N4+N+mem2edge,:);
    tempload(N+mem2edge) = tauelem(1:Nmem-1,1).*tempload(N+mem2edge)+tauelem(1:Nmem-1,2).*tempload(N4+N+mem2edge);             

    tempStiff(N4+N+mem2edge,:) = sparse([1:Nmem-1 1:Nmem-1]',[N+mem2edge; N4+N+mem2edge],[nelem(:,1); nelem(:,2)],Nmem-1,2*N4+N3+2*N2+N5);
    tempload(N4+N+mem2edge) = dydt(Nmem+1:2*Nmem-1).*nelem(:,2);

%% --------------------------------------------------------------------------------
    % extract unknown variables' indices
 
    % bdc of fluid
    bdnodeup = geom.bdnodeup; bdedgeup = geom.bdedgeup;

    ind1x = [bdnodeup;bdedgeup]; ind1x(geom.mem2node(floor((Nmem+1)/2))) = 1; 
    ind1y = [bdnodeup;bdedgeup];
    ind1p = zeros(N3,1); ind1p(1) = 1;

    ind3 = [0;zeros(Nffi-2,1);0;zeros(Nffi-1,1)]; ind2D2 = [1;zeros(Nmem-2,1);1;0;0;zeros(Nmem-1,1)];

    indexdelete = [ind1x;ind1y;ind1p;ind3;ind3;ind2D2];
    
    Stiffmat = tempStiff(indexdelete(:)==0, indexdelete(:)==0);
    
%     condest(Stiffmat)

    load = tempload(indexdelete(:)==0);
    
    solved = Stiffmat\load;
    
%% -----------------------------------------------------------------------------
    % arrange solutions
    back = zeros(size(indexdelete));
    back(indexdelete(:)==0)=solved;

    flu_vel = [back(1:N4) back(N4+1:2*N4)];
    sln.flu_vel = flu_vel;
    pressure = back(2*N4+1:2*N4+N3);
    sln.pressure = pressure;

    sln.pffi = pffi+[back(2*N4+N3+1:2*N4+N3+N2) back(2*N4+N3+N2+1:2*N4+N3+2*N2)];
    
    sln.nu = back(2*N4+N3+2*N2+1:2*N4+N3+2*N2+N5);
     
    % visualization
    visual = false;
    while visual
        node(geom.mem2node,2) = sln.pmem(1:para.Nmem,2);
        pos = [node; 0.5*(node(geom.edge(:,1),1:2)+node(geom.edge(:,2),1:2))];
        trimesh(elem,node(:,1),node(:,2),zeros(size(node,1),1));
        hold on;
        quiver(pos(:,1),pos(:,2),flu_vel(:,1),flu_vel(:,2),'r');
        plot(pmem(cl2mem1,1),pmem(cl2mem1,2),'rs');
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
%         disp(['nu symmetry check: ',num2str(max(abs(sln.nu(1:Nmem)-sln.nu(Nmem:-1:1))))]);
        %---------------------------------------------------------------
        
        plot(sln.pmem(1:Nmem,1),sln.pmem(1:Nmem,2),'o-r');
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
        
%         plot(sln.pmem(1:Nmem,1),sln.pmem(1:Nmem,2),'o-r');
%         hold on;
%         plot(sln.pffi(1:Nffi,1),sln.pffi(1:Nffi,2),'o-r');
%         plot(sln.pmemold(1:Nmem,1),sln.pmemold(1:Nmem,2),'*--b');
%         plot(sln.pffiold(1:Nffi,1),sln.pffiold(1:Nffi,2),'*--b');
%         axis equal;
%         hold off; 
        
        plot(pmem_nodal(1:cl2mem1,1),sln.nu(1:cl2mem1),'-ob');
        hold on;
        plot(pmem_nodal(cl2mem1:cl2mem2,1),[sln.nu(Nmem+1);sln.nu(cl2mem1+1:cl2mem2-1);sln.nu(Nmem+2)],'-sr');
        plot(pmem_nodal(cl2mem2:Nmem,1),sln.nu(cl2mem2:Nmem),'-ob');    
        hold off;     
%         pause(0.001);
        
%         disp(['left nu jump: ', num2str(sln.nu(Nmem+1)-sln.nu(cl2mem1))]);    
%         disp(['right nu jump: ', num2str(sln.nu(Nmem+2)-sln.nu(cl2mem2))]); 
   
        break;
    end
    
        plot(sln.pmem(1:Nmem,1),sln.pmem(1:Nmem,2),'o-r');
        hold on;
        plot(sln.pffi(1:Nffi,1),sln.pffi(1:Nffi,2),'o-r');
        plot(sln.pmemold(1:Nmem,1),sln.pmemold(1:Nmem,2),'*--b');
        plot(sln.pffiold(1:Nffi,1),sln.pffiold(1:Nffi,2),'*--b');
        axis equal;
        hold off;  
%         pause(0.01);
    
     sln.temp_pmem = sln.pmem;
     sln.pmem(cl2mem1,:) = sln.pffi(1,:);    sln.pmem(cl2mem2,:) = sln.pffi(Nffi,:);
    
     para.curlen = getlength(sln.pmem(1:Nmem,:));
     bd_droplet = [sln.pffi(1:Nffi,:); sln.pmem(cl2mem2-1:-1:cl2mem1+1,:)];
     areanow = polyarea(bd_droplet(:,1),bd_droplet(:,2));
     para.areanow = areanow;
     disp(['droplet area after ffi: ',num2str(areanow)]);
     
     disp(['symmetry check cl: ',num2str(sln.pmem(cl2mem1,1)+sln.pmem(cl2mem2,1))]);
    

end