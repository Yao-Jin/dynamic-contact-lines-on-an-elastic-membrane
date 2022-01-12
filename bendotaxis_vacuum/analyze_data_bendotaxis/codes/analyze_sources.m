clear;
% clc;
addpath(genpath('.'));
dbstop if error

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end

%1 4 9 17
% load('record_wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1_ind_17.mat'); 

%1 4 9 14
load('record_dewetting_20_10_gap_4_xl_2_gamma2_0_3_cb_5e3_Ca_0_1_ind_14.mat');

disp(['Time = ',num2str((hstry.itr-1)*para.dt)]);

node = geom.node; elem = geom.elem; vel = sln.flu_vel;
pressure = sln.pressure; N = length(node); pffil = sln.pffilold; pffir = sln.pffirold;
Nffil = para.Nffil; Nffir = para.Nffir;
cb = para.cb; gamma1 = para.gamma1; Ca = para.Ca;
nu = sln.nu; kappa = sln.kmem;
pmem = sln.pmemold; Nmem = para.Nmem; Nm = para.Nm; x = pmem(:,1); y = pmem(:,2);
cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2);
[nelem, nnode, tauelem, taunode] = getnormal_tau_mem_P1(pmem(1:Nmem,1:2));
s_elem = sqrt((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2);
pypxelem = (pmem(2:Nmem,2)-pmem(1:Nmem-1,2))./(pmem(2:Nmem,1)-pmem(1:Nmem-1,1));
lenelem = sqrt(1+pypxelem.^2);
pypxnode = -nnode(:,1)./nnode(:,2);
lenelem_node = sqrt(1+pypxnode.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gap = 10;
% trimesh(elem,node(:,1),node(:,2)-para.By,pressure(1:N),'facecolor','interp','edgecolor','interp'); 
% hold on;
% trimesh(elem,node(:,1),para.By-node(:,2),pressure(1:N),'facecolor','interp','edgecolor','interp');
% plot(pffil(1:Nffil,1),pffil(1:Nffil,2)-para.By,'-r','LineWidth',1.0);
% plot(pffil(1:Nffil,1),para.By-pffil(1:Nffil,2),'-r','LineWidth',1.0);
% plot(pffir(1:Nffir,1),pffir(1:Nffir,2)-para.By,'-r','LineWidth',1.0);
% plot(pffir(1:Nffir,1),para.By-pffir(1:Nffir,2),'-r','LineWidth',1.0);
% plot(pmem(cl2mem1-gap:cl2mem2+gap,1),pmem(cl2mem1-gap:cl2mem2+gap,2)-para.By,'-r','LineWidth',1.0);
% plot(pmem(cl2mem1-gap:cl2mem2+gap,1),para.By-pmem(cl2mem1-gap:cl2mem2+gap,2),'-r','LineWidth',1.0);
% plot([0 0],[-para.By, para.By],'-k','LineWidth',3.0); 
% 
% quiver(node(:,1),node(:,2)-para.By,vel(1:N,1),vel(1:N,2),'k','LineWidth',1.5,'AutoScaleFactor',0.6);
% quiver(node(:,1),para.By-node(:,2),vel(1:N,1),-vel(1:N,2),'k','LineWidth',1.5,'AutoScaleFactor',0.6);
% 
% view(2);
% axis equal;
% xlim([pmem(cl2mem1-gap,1), pmem(cl2mem2+gap,1)]); ylim([-para.By-0.05, para.By+0.05]);
% colorbar('northoutside');
% set(gca,'FontSize',32);
% set(gca,'YAxisLocation','right');
% hold off;

% min(vel(1:N,1))
% sum(vel(1:N,1))/N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot bragamma
% plot(x(cl2mem1:cl2mem2),nu(1:Nm),'-ob');
% hold on;
% plot([x(1) x(cl2mem1)],[nu(Nm+1) nu(Nm+1)],'-r');
% plot([x(cl2mem2) x(Nmem)],[0 0],'-r');
% hold off;

% plot dnuds
dnuds = zeros(Nm,1);
dnuds(1) = (-3*nu(1)+4*nu(2)-nu(3))/(2*s_elem(cl2mem1));
dnuds(2:Nm-1) = (nu(3:Nm)-nu(1:Nm-2))./(s_elem(cl2mem1:cl2mem2-2)+s_elem(cl2mem1+1:cl2mem2-1));
dnuds(Nm) = (3*nu(Nm)-4*nu(Nm-1)+nu(Nm-2))/(2*s_elem(cl2mem2-1));

% kappa = getcurvature(pmem(1:Nmem,1:2));

% plot(x(cl2mem1-1:cl2mem2+1), pypxnode(cl2mem1-1:cl2mem2+1), '-or');

% compute dkappads
dkappads = zeros(Nmem+2,1);
dkappads(1) = (-3*kappa(1)+4*kappa(2)-kappa(3))/(2*s_elem(1));
dkappads(2:cl2mem1-1) = (kappa(3:cl2mem1)-kappa(1:cl2mem1-2))./(s_elem(1:cl2mem1-2)+s_elem(2:cl2mem1-1));
dkappads(cl2mem1) = (3*kappa(cl2mem1)-4*kappa(cl2mem1-1)+kappa(cl2mem1-2))/(2*s_elem(cl2mem1-1));

dkappads(Nmem+1) = (-3*kappa(cl2mem1)+4*kappa(cl2mem1+1)-kappa(cl2mem1+2))/(2*s_elem(cl2mem1));
dkappads(cl2mem1+1:cl2mem2-1) = (kappa(cl2mem1+2:cl2mem2)-kappa(cl2mem1:cl2mem2-2))./(s_elem(cl2mem1:cl2mem2-2)+s_elem(cl2mem1+1:cl2mem2-1));
dkappads(Nmem+2) = (3*kappa(cl2mem2)-4*kappa(cl2mem2-1)+kappa(cl2mem2-2))/(2*s_elem(cl2mem2-1));

dkappads(cl2mem2) = (-3*kappa(cl2mem2)+4*kappa(cl2mem2+1)-kappa(cl2mem2+2))/(2*s_elem(cl2mem2));
dkappads(cl2mem2+1:Nmem-1) = (kappa(cl2mem2+2:Nmem)-kappa(cl2mem2:Nmem-2))./(s_elem(cl2mem2:Nmem-2)+s_elem(cl2mem2+1:Nmem-1));
dkappads(Nmem) = (3*kappa(Nmem)-4*kappa(Nmem-1)+kappa(Nmem-2))/(2*s_elem(Nmem-1));

% compute d2kappads
d2kappads = zeros(Nmem,1);


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


source1 = cb*d2kappads.*nnode(cl2mem1:cl2mem2,1);

source2 = cb/2*kappa(cl2mem1:cl2mem2).^3.*nnode(cl2mem1:cl2mem2,1);

source3 = (nu(1:Nm)-gamma1).*kappa(cl2mem1:cl2mem2).*nnode(cl2mem1:cl2mem2,1);

source4 = dnuds(1:Nm).*taunode(cl2mem1:cl2mem2,1);

% plot(x(cl2mem1:cl2mem2),source1,'-or');
% hold on;
% plot(x(cl2mem1:cl2mem2),source4,'-sb');
% hold off;

costhetadL = nu(1)-nu(Nm+1)-para.gamma1+para.gamma2;
costhetadR = nu(Nm)-para.gamma1+para.gamma2;
% sinthetadL = cb*(dkappads(Nmem+1)-dkappads(cl2mem1))
% sinthetadR = -cb*(dkappads(Nmem+2)-dkappads(cl2mem2))
sinthetadL = sqrt(1-costhetadL^2);
sinthetadR = sqrt(1-costhetadR^2);

    f_nor = 1/Ca*(source1+source2);
    f_tan = -1/Ca*(source4);

    F_nor = sum(1/2*(f_nor(1:Nm-1)+f_nor(2:Nm)).*s_elem(cl2mem1:cl2mem2-1))
    F_tan = sum(1/2*(f_tan(1:Nm-1)+f_tan(2:Nm)).*s_elem(cl2mem1:cl2mem2-1));
%     Fcl_l = 1/Ca*(sin(para.thetaY)*nnode(cl2mem1,1)-cos(para.thetaY)*taunode(cl2mem1,1))
%     Fcl_r = 1/Ca*(sin(para.thetaY)*nnode(cl2mem2,1)+cos(para.thetaY)*taunode(cl2mem2,1))
    Fcl_l = 1/Ca*(sinthetadL*nnode(cl2mem1,1)-costhetadL*taunode(cl2mem1,1));
    Fcl_r = 1/Ca*(sinthetadR*nnode(cl2mem2,1)+costhetadR*taunode(cl2mem2,1));
    Fcl = Fcl_l+Fcl_r
%     horizon_join_force = F_nor+F_tan+Fcl_l+Fcl_r


% horizon = source1-source4+source2+source3;
% 
% % numerical integration
% integral = sum(1/2*(horizon(1:Nm-1)+horizon(2:Nm)).*s_elem(cl2mem1:cl2mem2-1))
% % getlength(pmem(cl2mem1:cl2mem2,1:2))
% 
% % outnorl = pffil(1,1:2)-pffil(2,1:2); outnorr = pffir(1,1:2)-pffir(2,1:2);
% % outnorl = outnorl./sqrt(outnorl(1)^2+outnorl(2)^2);
% % outnorr = outnorr./sqrt(outnorr(1)^2+outnorr(2)^2);
% % thetadl = acos(-outnorl*taunode(cl2mem1,1:2)')
% % thetadr = acos(outnorr*taunode(cl2mem2,1:2)')
% % sin(thetadl)*nnode(cl2mem1,1)+sin(thetadr)*nnode(cl2mem2,1)
% % source_cl = para.gamma3*(outnorr(1)+outnorl(1))
% 
% % thetaY = para.thetaY
% 
% source_cl = sin(para.thetaY)*nnode(cl2mem1,1)+sin(para.thetaY)*nnode(cl2mem2,1)
% 
% source5 = -cos(para.thetaY)*taunode(cl2mem1,1)+cos(para.thetaY)*taunode(cl2mem2,1)
% 
% horizonal_join_force = integral+source_cl+source5
% 
% horizonal_join_force/integral
% 
% % thetad = acos(para.gamma2-para.gamma1+nu(1)-nu(Nm+1))
% 
% % sin(thetad)*nnode(cl2mem1,1)+sin(thetad)*nnode(cl2mem2,1)
% 
% D_dkappads = (-3*kappa(cl2mem1)+4*kappa(cl2mem1+1)-kappa(cl2mem1+2))/(2*s_elem(cl2mem1))...
%              -(3*kappa(cl2mem1)-4*kappa(cl2mem1-1)+kappa(cl2mem1-2))/(2*s_elem(cl2mem1-1));         
% D_dkappads = cb*D_dkappads;
% 
% % thetad = asin(D_dkappads)
% 
% % sin(thetad)*nnode(cl2mem1,1)+sin(thetad)*nnode(cl2mem2,1)
