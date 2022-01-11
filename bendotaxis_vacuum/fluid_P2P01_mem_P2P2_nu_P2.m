function [sln, para] = fluid_P2P01_mem_P2P2_nu_P2(geom,para,sln)

%  construct and solve the coupled system

    elem = geom.elem; node = geom.node; edge = geom.edge;
    Nt = size(elem,1); N = size(node,1); Ne = size(edge,1);
    cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2);
    Nmem = para.Nmem; Nffil = para.Nffil; Nffir = para.Nffir; dt = para.dt;
    pmem = sln.pmem(:,1:2); pffil = sln.pffil(:,1:2); pffir = sln.pffir(:,1:2); 
    
% construct the stiffness matrix and the load vector by assembling local stiffness matrices

% stiffness matrix for bulk elements
    [stiff_bulk] = stiff_Bulk_P2P01(geom,para,sln);
    B11 = stiff_bulk.B11; B12 = stiff_bulk.B12; B13 = stiff_bulk.B13; B14 = stiff_bulk.B14;
    B21 = stiff_bulk.B21; B22 = stiff_bulk.B22;
    
% stiffness matrix for fluid interface elements
    [stiff_interface] = stiff_Interface_P2(geom,para,sln);
    Il1 = stiff_interface.Il1; Ir1 = stiff_interface.Ir1;

% stiffness matrix for membrane elements
    [stiff_membrane] = stiff_Membrane_P2P2_vacuum(geom,para,sln);
    M11 = stiff_membrane.M11; M12 = stiff_membrane.M12; M13 = stiff_membrane.M13; M14 = stiff_membrane.M14;
    M15 = stiff_membrane.M15; M16 = stiff_membrane.M16; M17 = stiff_membrane.M17; M18 = stiff_membrane.M18; M19 = stiff_membrane.M19; 
    M21 = stiff_membrane.M21; M22 = stiff_membrane.M22; M23 = stiff_membrane.M23;
    M31 = stiff_membrane.M31; M32 = stiff_membrane.M32; M33 = stiff_membrane.M33;
    M41 = stiff_membrane.M41; M42 = stiff_membrane.M42;
 
    [stiff_nu] = stiff_nu_P2_vacuum(geom,para,sln);
    G11 = stiff_nu.G11; G12 = stiff_nu.G12; 
    G13 = stiff_nu.G13; G14 = stiff_nu.G14; G15 = stiff_nu.G15; 
    G2 = stiff_nu.G2;
    G31 = stiff_nu.G31; G32 = stiff_nu.G32;
    
% stiffness matrix for elements at the contact line  
    [stiff_cl] = stiff_Contactline_P2(geom,para,sln);
    cl1 = stiff_cl.cl1; cl2 = stiff_cl.cl2; cl3 = stiff_cl.cl3; cl4 = stiff_cl.cl4;
    cl7 = stiff_cl.cl7; cl8 = stiff_cl.cl8;
    cl9 = stiff_cl.cl9; cl10 = stiff_cl.cl10; cl11 = stiff_cl.cl11;
    cl12 = stiff_cl.cl12; cl13 = stiff_cl.cl13; 
    CL1 = stiff_cl.CL1;
    % u_x u_y p y kappa nu
    % here we are solving the difference of y from time t^m to t^(m+1)
    Nm = cl2mem2-cl2mem1+1;
    N1 = 2*Nmem-1; N3 = Nt+N; N4 = N+Ne; N5 = 2*Nm;  
             
    tempStiff = [B11 B12 -B21' -M11 M13 -G11-cl3-G31;
                 B13 B14 -B22' -M12 M14 -G12-cl4-G32;
                 B21 B22 sparse(N3,N3+2*N1+N5);
                 M18 M19 sparse(N1,N3) M17/dt-M21-M31 M22+M32 -G2;
                 sparse(N1,2*N4+N3) M42 M41 sparse(N1,N5);
                 G14+cl9 G15+cl10 sparse(N5,N3) -cl8/dt sparse(N5,N1) G13+cl7+cl11];             
             
%     tempload = [-cl1-cl12;-cl2-cl13;sparse(N3,1);sparse(N1,1);sparse(N1+N5,1)];
    tempload = [-CL1;sparse(N4,1);sparse(N3,1);sparse(N1,1);sparse(N1+N5,1)];
    
    tempload(1:N4) = tempload(1:N4)-Il1*pffil(:,1)-Ir1*pffir(:,1)+(M11-M15)*pmem(:,2);
    tempload(N4+1:2*N4) = tempload(N4+1:2*N4)-Il1*pffil(:,2)-Ir1*pffir(:,2)+(M12-M16)*pmem(:,2);
    
    tempload(2*N4+N3+1:2*N4+N3+N1) = tempload(2*N4+N3+1:2*N4+N3+N1)+(M21+M31-M23-M33)*pmem(:,2);
    
    tempload(2*N4+N3+N1+1:2*N4+N3+2*N1) = tempload(2*N4+N3+N1+1:2*N4+N3+2*N1)-M42*pmem(:,2);     
         
%% --------------------------------------------------------------------------------
    % extract unknown variables' indices
 
    % bdc of fluid
    bdnodeup = geom.bdnodeup; bdedgeup = geom.bdedgeup;
    
    ind1x = zeros(N4,1);
    ind1y = [bdnodeup;bdedgeup];
    ind1p = zeros(N3,1); ind1p(Nt+1) = 1;
    ind2L = [1;zeros(Nmem-2,1);0;zeros(Nmem-1,1)]; ind2R = [0;zeros(Nmem-2,1);1;zeros(Nmem-1,1)];
    ind3 = zeros(2*Nm,1);

    indexdelete = [ind1x;ind1y;ind1p;ind2L;ind2R;ind3];

    Stiffmat = tempStiff(indexdelete(:)==0, indexdelete(:)==0);
    
%     condest(Stiffmat)
    
    load = tempload(indexdelete(:)==0);
    
    solved = Stiffmat\load;
    
%% -----------------------------------------------------------------------------
    % arrange solutions
    back = zeros(size(indexdelete));
    back(indexdelete(:)==0)=solved;

    sln.pmemold = sln.pmem; sln.kmemold = sln.kmem; sln.nu_old = sln.nu; 
    sln.pffilold = sln.pffil; sln.pffirold = sln.pffir;
    
    flu_vel = [back(1:N4) back(N4+1:2*N4)];
    sln.flu_vel = flu_vel;
    pressure = back(2*N4+1:2*N4+N3);
    sln.pressure = pressure;
    
    index = 2*N4+N3;
    sln.pmem(:,2) = pmem(:,2)+back(index+1:index+N1); sln.maxdy = max(abs(back(index+1:index+N1)));
    sln.dydt = back(index+1:index+N1)/para.dt;
    
    sln.kmem = back(index+N1+1:index+2*N1);
    sln.nu = back(index+2*N1+1:index+2*N1+N5);    
    
    sln.pffil(1,2) = sln.pmem(cl2mem1,2); sln.pffir(1,2) = sln.pmem(cl2mem2,2);
      
    % visualization
    visual = true;
    while visual
        pos = [node; 0.5*(node(geom.edge(:,1),1:2)+node(geom.edge(:,2),1:2))];
        trimesh(elem,node(:,1),node(:,2),zeros(size(node,1),1));
        hold on;
        quiver(pos(:,1),pos(:,2),flu_vel(:,1),flu_vel(:,2),'r');
        plot(pmem(cl2mem1,1),pmem(cl2mem1,2),'bo');
        view(2);
        axis equal;
        hold off;
%         pause(0.01);
               
%         elem_mid(:,1) = 1/3*(node(elem(:,1),1)+node(elem(:,2),1)+node(elem(:,3),1));
%         elem_mid(:,2) = 1/3*(node(elem(:,1),2)+node(elem(:,2),2)+node(elem(:,3),2));
%         plot3(elem_mid(:,1),elem_mid(:,2),pressure(N+1:N+Nt),'or'); 
        
        pmem = sln.pmem;  pffil = sln.pffil; pffir = sln.pffir;
        plot(pmem(1:Nmem,1),pmem(1:Nmem,2),'o-r');
        hold on;
        plot(pffil(1:Nffil,1),pffil(1:Nffil,2),'o-r');
        plot(pffir(1:Nffir,1),pffir(1:Nffir,2),'o-r');
        plot(sln.pmemold(1:Nmem,1),sln.pmemold(1:Nmem,2),'*--b');
        plot(sln.pffilold(1:Nffil,1),sln.pffilold(1:Nffil,2),'*--b');
        plot(sln.pffirold(1:Nffir,1),sln.pffirold(1:Nffir,2),'*--b');
        axis equal;
        hold off;
        
        plot(pmem(1:Nmem,1),sln.kmem(1:Nmem),'-*r');
        
        pmem = sln.pmemold(1:Nmem,1:2);
        plot(pmem(cl2mem1:cl2mem2,1),sln.nu(1:Nm),'-sr');
        hold on;
        plot([pmem(1,1) pmem(cl2mem1,1)],[sln.nu(Nm+1) sln.nu(Nm+1)],'-ob');
        plot([pmem(cl2mem2,1) pmem(Nmem,1)],[0 0],'-ob');
        hold off;     
%         pause(0.001);
        
%         disp(['left nu jump: ', num2str(sln.nu(1)-sln.nu(Nm+1))]);    
%         disp(['right nu jump: ', num2str(sln.nu(Nm))]);      

        break;
    end
     
     para.curlen = getlength(sln.pmem(1:Nmem,:));
     bd_droplet = [sln.pffil(1:Nffil,1:2); sln.pffir(Nffir:-1:1,1:2); sln.pmem(cl2mem2-1:-1:cl2mem1+1,1:2)];
     areanow = polyarea(bd_droplet(:,1),bd_droplet(:,2));     para.areanow = areanow;
     disp(['droplet area after mem: ',num2str(areanow)]);     

end