function [sln, para] = fluid_P2P01_ffi_P2_nu_P2(geom,para,sln)

%  construct and solve the coupled system

    elem = geom.elem; node = geom.node; edge = geom.edge;
    Nt = size(elem,1); N = size(node,1); Ne = size(edge,1);
    cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2);
    Nmem = para.Nmem; Nffil = para.Nffil; Nffir = para.Nffir; Nm = para.Nm;
    pmem = sln.pmem; pffil = sln.pffil; pffir = sln.pffir; dt = para.dt;
    
% construct the stiffness matrix and the load vector by assembling local stiffness matrices

% stiffness matrix for bulk elements
    [stiff_bulk] = stiff_Bulk_P2P01(geom,para,sln);
    B11 = stiff_bulk.B11; B12 = stiff_bulk.B12; B13 = stiff_bulk.B13; B14 = stiff_bulk.B14;
    B21 = stiff_bulk.B21; B22 = stiff_bulk.B22;
    
% stiffness matrix for fluid interface elements
    [stiff_interface] = stiff_Interface_P2(geom,para,sln);
    Il1 = stiff_interface.Il1; Il2 = stiff_interface.Il2; Il3 = stiff_interface.Il3;
    Ir1 = stiff_interface.Ir1; Ir2 = stiff_interface.Ir2; Ir3 = stiff_interface.Ir3;
    
    [stiff_nu] = stiff_nu_P2_vacuum(geom,para,sln);
    G11 = stiff_nu.G11; G12 = stiff_nu.G12; 
    G13 = stiff_nu.G13; G14 = stiff_nu.G14; G15 = stiff_nu.G15; 
    
% stiffness matrix for elements at the contact line  
    [stiff_cl] = stiff_Contactline_P2(geom,para,sln);
    cl1 = stiff_cl.cl1; cl2 = stiff_cl.cl2; cl3 = stiff_cl.cl3; cl4 = stiff_cl.cl4;
    cl7 = stiff_cl.cl7; cl8 = stiff_cl.cl8;
    cl9 = stiff_cl.cl9; cl10 = stiff_cl.cl10; cl11 = stiff_cl.cl11;
    
    % u_x u_y p r_x r_y nu
    % here we are solving the difference of r_x, r_y from time t^m to t^(m+1)
    N2l = 2*Nffil-1; N2r = 2*Nffir-1; N3 = Nt+N; N4 = N+Ne; N5 = 2*Nm;
    
    tempStiff = [B11 B12 -B21' Il1 sparse(N4,N2l) Ir1 sparse(N4,N2r) -G11-cl3;
                 B13 B14 -B22' sparse(N4,N2l) Il1 sparse(N4,N2r) Ir1 -G12-cl4;
                 B21 B22 sparse(N3,N3+2*N2l+2*N2r+N5);
                 -Il3 sparse(N2l,N4+N3) Il2/dt sparse(N2l,N2l+2*N2r+N5);
                 sparse(N2l,N4) -Il3 sparse(N2l,N3+N2l) Il2/dt sparse(N2l,2*N2r+N5);
                 -Ir3 sparse(N2r,N4+N3+2*N2l) Ir2/dt sparse(N2r,N2r+N5);
                 sparse(N2r,N4) -Ir3 sparse(N2r,N3+2*N2l+N2r) Ir2/dt sparse(N2r,N5);
                 G14+cl9 G15+cl10 sparse(N5,N3+2*N2l+2*N2r) G13+cl7+cl11];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
    Dy1 = 1./(pmem(cl2mem1+1,1)-pmem(cl2mem1-1,1))*[pmem(cl2mem1+1,1)-pmem(cl2mem1,1) pmem(cl2mem1,1)-pmem(cl2mem1-1,1)]*...
          [(pmem(cl2mem1,2)-pmem(cl2mem1-1,2))/(pmem(cl2mem1,1)-pmem(cl2mem1-1,1));(pmem(cl2mem1+1,2)-pmem(cl2mem1,2))/(pmem(cl2mem1+1,1)-pmem(cl2mem1,1))];

    Dy2 = 1./(pmem(cl2mem2+1,1)-pmem(cl2mem2-1,1))*[pmem(cl2mem2+1,1)-pmem(cl2mem2,1) pmem(cl2mem2,1)-pmem(cl2mem2-1,1)]*...
          [(pmem(cl2mem2,2)-pmem(cl2mem2-1,2))/(pmem(cl2mem2,1)-pmem(cl2mem2-1,1));(pmem(cl2mem2+1,2)-pmem(cl2mem2,2))/(pmem(cl2mem2+1,1)-pmem(cl2mem2,1))];
      
    tempStiff(2*N4+N3+N2l+1,:) = [sparse(1,2*N4+N3) -Dy1/dt sparse(1,N2l-1) 1/dt sparse(1,N2l-1) sparse(1,2*N2r+N5)];
    tempStiff(2*N4+N3+2*N2l+N2r+1,:) = [sparse(1,2*N4+N3+2*N2l) -Dy2/dt sparse(1,N2r-1) 1/dt sparse(1,N2r-1+N5)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    tempload = [-cl1;-cl2;sparse(N3+2*N2l+2*N2r+N5,1)];
                
    tempload(1:N4) = tempload(1:N4)-Il1*pffil(:,1)-Ir1*pffir(:,1);
    tempload(N4+1:2*N4) = tempload(N4+1:2*N4)-Il1*pffil(:,2)-Ir1*pffir(:,2);
                
    tmp = 2*N4+N3+2*N2l+2*N2r;
    tempload(tmp+1:tmp+N5) = tempload(tmp+1:tmp+N5)+cl8*sln.dydt;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    [nelem, nnode, tauelem, taunode] = getnormal_tau_mem_P1(pmem(1:Nmem,1:2));
    dydt = sln.dydt; mem2node = geom.mem2node; mem2edge = geom.mem2edge;
    
    index1 = (cl2mem1:cl2mem2);
    tempStiff(mem2node,:) = taunode(index1,1).*tempStiff(mem2node,:)...
                          + taunode(index1,2).*tempStiff(N4+mem2node,:);
    tempload(mem2node) = taunode(index1,1).*tempload(mem2node)+taunode(index1,2).*tempload(N4+mem2node);             
    
    index2 = (1:cl2mem2-cl2mem1+1);
    tempStiff(N4+mem2node,:) = sparse([index2 index2]',[mem2node; N4+mem2node]', [nnode(index1,1); nnode(index1,2)],...
                                       cl2mem2-cl2mem1+1,2*N4+N3+2*N2l+2*N2r+N5);
    tempload(N4+mem2node) = dydt(index1).*nnode(index1,2);    
    
    index1 = (cl2mem1:cl2mem2-1);
    tempStiff(N+mem2edge,:) = tauelem(index1,1).*tempStiff(N+mem2edge,:)...
                               + tauelem(index1,2).*tempStiff(N4+N+mem2edge,:);
    tempload(N+mem2edge) = tauelem(index1,1).*tempload(N+mem2edge)+tauelem(index1,2).*tempload(N4+N+mem2edge);             
    
    index2 = (1:cl2mem2-cl2mem1);
    tempStiff(N4+N+mem2edge,:) = sparse([index2 index2]',[N+mem2edge; N4+N+mem2edge],[nelem(index1,1); nelem(index1,2)],...
                                        cl2mem2-cl2mem1,2*N4+N3+2*N2l+2*N2r+N5);
    tempload(N4+N+mem2edge) = dydt(Nmem+cl2mem1:Nmem+cl2mem2-1).*nelem(index1,2);
    
%% --------------------------------------------------------------------------------
    % extract unknown variables' indices
 
    % bdc of fluid
    bdnodeup = geom.bdnodeup; bdedgeup = geom.bdedgeup;
    
    ind1x = zeros(N4,1);
    ind1y = [bdnodeup;bdedgeup];
    ind1p = zeros(N3,1); ind1p(Nt+1) = 1;
    ind3lx = [0;zeros(Nffil-2,1);0;zeros(Nffil-1,1)]; 
    ind3rx = [0;zeros(Nffir-2,1);0;zeros(Nffir-1,1)];
    ind3ly = [0;zeros(Nffil-2,1);1;zeros(Nffil-1,1)]; 
    ind3ry = [0;zeros(Nffir-2,1);1;zeros(Nffir-1,1)];
    ind3 = zeros(2*Nm,1);
    indexdelete = [ind1x;ind1y;ind1p;ind3lx;ind3ly;ind3rx;ind3ry;ind3];
    
    Stiffmat = tempStiff(indexdelete(:)==0, indexdelete(:)==0);
    
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

    sln.pffil = pffil+[back(2*N4+N3+1:2*N4+N3+N2l) back(2*N4+N3+N2l+1:2*N4+N3+2*N2l)];
    sln.pffir = pffir+[back(2*N4+N3+2*N2l+1:2*N4+N3+2*N2l+N2r) back(2*N4+N3+2*N2l+N2r+1:2*N4+N3+2*N2l+2*N2r)];
    sln.nu = back(2*N4+N3+2*N2l+2*N2r+1:2*N4+N3+2*N2l+2*N2r+N5); 
    
    % visualization
    visual = true;
    while visual
        node(geom.mem2node,2) = sln.pmem(cl2mem1:cl2mem2,2);
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
        plot(pmem(:,1),pmem(:,2),'o-r');
        hold on;
        plot(pffil(:,1),pffil(:,2),'o-r');
        plot(pffir(:,1),pffir(:,2),'o-r');
        plot(sln.pmemold(1:Nmem,1),sln.pmemold(1:Nmem,2),'*--b');
        plot(sln.pffilold(1:Nffil,1),sln.pffilold(1:Nffil,2),'*--b');
        plot(sln.pffirold(1:Nffir,1),sln.pffirold(1:Nffir,2),'*--b');
        axis equal;
        hold off;  
%         pause(0.001);
        
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
%         pause(0.01);    
        

     sln.pmem(cl2mem1,:) = sln.pffil(1,:); sln.pmem(cl2mem2,:) = sln.pffir(1,:);
     para.curlen = getlength(sln.pmem(1:Nmem,:));
     bd_droplet = [sln.pffil(1:Nffil,1:2); sln.pffir(Nffir:-1:1,1:2); sln.pmem(cl2mem2-1:-1:cl2mem1+1,1:2)];
     areanow = polyarea(bd_droplet(:,1),bd_droplet(:,2));     para.areanow = areanow;
     disp(['droplet area after ffi: ',num2str(areanow)]);     

end