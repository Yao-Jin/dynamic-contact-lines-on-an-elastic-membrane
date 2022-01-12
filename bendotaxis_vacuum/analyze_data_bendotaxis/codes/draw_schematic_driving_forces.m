clear;
clc;
addpath(genpath('.'));
dbstop if error

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end

load('record_wetting_gap_1_xl_1_gamma1_0_5_cb_5e1_Ca_1_t_8.mat');
% load('record_dewetting_gap_2_xl_0_5_gamma2_0_1_cb_1e1_Ca_1_t_14.mat');

pffil = sln.pffilold; pffir = sln.pffirold; Nffil = para.Nffil; Nffir = para.Nffir;
cb = para.cb; gamma1 = para.gamma1; Ca = para.Ca; By = para.By; Bx = para.Bx; thetaY = para.thetaY;
nu = sln.nu; kappa = sln.kmem;
pmem = sln.pmemold; Nmem = para.Nmem; Nm = para.Nm; x = pmem(:,1); y = pmem(:,2);
cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2);
[nelem, nnode, tauelem, taunode] = getnormal_tau_mem_P1(pmem(1:Nmem,1:2));
s_elem = sqrt((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2);
pypxelem = (pmem(2:Nmem,2)-pmem(1:Nmem-1,2))./(pmem(2:Nmem,1)-pmem(1:Nmem-1,1));
lenelem = sqrt(1+pypxelem.^2);
pypxnode = -nnode(:,1)./nnode(:,2);
lenelem_node = sqrt(1+pypxnode.^2);

% compute dnuds
thetadL = acos(nu(1)-nu(Nm+1)+para.gamma2-para.gamma1);
thetadR = acos(nu(Nm)-0+para.gamma2-para.gamma1);
dnuds = zeros(Nm,1);
dnuds(1) = (-3*nu(1)+4*nu(2)-nu(3))/(2*s_elem(cl2mem1));
dnuds(2:Nm-1) = (nu(3:Nm)-nu(1:Nm-2))./(s_elem(cl2mem1:cl2mem2-2)+s_elem(cl2mem1+1:cl2mem2-1));
dnuds(Nm) = (3*nu(Nm)-4*nu(Nm-1)+nu(Nm-2))/(2*s_elem(cl2mem2-1));

% compute d2kappads
d2kappads = zeros(Nm,1);
d2kappads(1) = ...
    (2*kappa(cl2mem1)-5*kappa(cl2mem1+1)+4*kappa(cl2mem1+2)-kappa(cl2mem1+3))...
    ./(s_elem(cl2mem1)^2);
d2kappads(2:Nm-1) = ...
    (kappa(cl2mem1:cl2mem2-2)+kappa(cl2mem1+2:cl2mem2)-2*kappa(cl2mem1+1:cl2mem2-1))...
    ./(((s_elem(cl2mem1:cl2mem2-2)+s_elem(cl2mem1+1:cl2mem2-1))/2).^2);
d2kappads(Nm) = ...
    (2*kappa(cl2mem2)-5*kappa(cl2mem2-1)+4*kappa(cl2mem2-2)-kappa(cl2mem2-3))...
    ./(s_elem(cl2mem2-1)^2);

source1 = cb*d2kappads;
source2 = cb/2*kappa(cl2mem1:cl2mem2).^3;
source3 = (nu(1:Nm)-gamma1).*kappa(cl2mem1:cl2mem2);
source4 = dnuds(1:Nm);
% normal_force = source1+source2+source3;
bend_force = source1+source2;

% numerical integration
% join_normal = sum(1/2*(normal_force(1:Nm-1)+normal_force(2:Nm)).*s_elem(cl2mem1:cl2mem2-1));
% join_tangent = sum(1/2*(source4(1:Nm-1)+source4(2:Nm)).*s_elem(cl2mem1:cl2mem2-1));

source_cl = sin(para.thetaY)*nnode(cl2mem1,1)+sin(para.thetaY)*nnode(cl2mem2,1);
source5 = -cos(para.thetaY)*taunode(cl2mem1,1)+cos(para.thetaY)*taunode(cl2mem2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw = 1.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start = 1;
plot(x(start:Nmem),y(start:Nmem)-By,'-','LineWidth',lw,'color',c2);
hold on;
plot(pffil(1:Nffil,1),pffil(1:Nffil,2)-By,'-','LineWidth',lw,'color',c1);
plot(pffir(1:Nffir,1),pffir(1:Nffir,2)-By,'-','LineWidth',lw,'color',c1);

% plot(x(start:Nmem),By-y(start:Nmem),'-','LineWidth',lw,'color',c2);
% plot(pffil(1:Nffil,1),By-pffil(1:Nffil,2),'-','LineWidth',lw,'color',c1);
% plot(pffir(1:Nffir,1),By-pffir(1:Nffir,2),'-','LineWidth',lw,'color',c1);

ind = cl2mem1+round((cl2mem2-cl2mem1)/2);
ind2 = cl2mem1+15; ind3 = cl2mem2-15;
pos = [x(ind) y(ind); x(cl2mem1) y(cl2mem1); x(cl2mem2) y(cl2mem2);
       x(ind2) y(ind2); x(ind3) y(ind3);
       x(ind) y(ind); x(ind2) y(ind2); x(ind3) y(ind3);
       x(cl2mem1) y(cl2mem1);
       x(cl2mem2) y(cl2mem2)];
   
% vector = [join_normal*nnode(ind,1:2); -join_tangent*taunode(ind,1:2)];
vector = [bend_force(ind-cl2mem1+1)*nnode(ind,1:2);
          sin(thetadL)*nnode(cl2mem1,1:2)-cos(thetadL)*taunode(cl2mem1,1:2);
          sin(thetadR)*nnode(cl2mem2,1:2)+cos(thetadR)*taunode(cl2mem2,1:2);
          bend_force(ind2-cl2mem1+1)*nnode(ind2,1:2);
          bend_force(ind3-cl2mem1+1)*nnode(ind3,1:2);
          bend_force(ind-cl2mem1+1)*nnode(ind,1) 0;
          bend_force(ind2-cl2mem1+1)*nnode(ind2,1) 0;
          bend_force(ind3-cl2mem1+1)*nnode(ind3,1) 0;
          sin(thetadL)*nnode(cl2mem1,1)-cos(thetadL)*taunode(cl2mem1,1) 0;
          sin(thetadR)*nnode(cl2mem2,1)+cos(thetadR)*taunode(cl2mem2,1) 0];

quiver(pos(:,1),pos(:,2)-By,vector(:,1),vector(:,2),'k','LineWidth',2,...
       'AutoScaleFactor',2,'MaxHeadSize',0.2);

% quiver(pos(:,1),By-pos(:,2),vector(:,1),-vector(:,2),'k','LineWidth',2,...
%        'AutoScaleFactor',0.8,'MaxHeadSize',0.1);
   
plot([pmem(start,1) pmem(start,1)],[pmem(start,2)-By 0],'-k','LineWidth',3.0);
% plot([pmem(start,1) pmem(start,1)],[pmem(start,2)+By 0],'-k','LineWidth',3.0);
plot([pmem(start,1),Bx],[0 0],'--k','LineWidth',1.0);
axis equal;
% xlim([0.7,2.8]);
% ylim([-1,1]);
axis off;
hold off;