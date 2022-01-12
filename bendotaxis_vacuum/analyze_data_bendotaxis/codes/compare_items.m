clear;
clc;
addpath(genpath('.'));
dbstop if error

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end

%1 4 9 17
load('record_wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1_ind_1.mat'); 

%1 4 9 14
% load('record_dewetting_20_10_gap_4_xl_2_gamma2_0_3_cb_5e3_Ca_0_1_ind_9.mat');

disp(['Time = ',num2str((hstry.itr-1)*para.dt)]);

cb1 = para.cb; gamma11 = para.gamma1; Ca1 = para.Ca;
nu1 = sln.nu; kappa1 = sln.kmem;
pmem1 = sln.pmemold; Nmem1 = para.Nmem; Nm1 = para.Nm; x1 = pmem1(1:Nmem1,1); y1 = pmem1(1:Nmem1,2);
cl2mem11 = geom.cl2mem(1); cl2mem12 = geom.cl2mem(2);
pressure1 = sln.pressure; node1 = geom.node; mem2node1 = geom.mem2node;
s_elem1 = sqrt((x1(2:end)-x1(1:end-1)).^2+(y1(2:end)-y1(1:end-1)).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute dnuds
dnuds1 = zeros(Nm1,1);
dnuds1(1) = (-3*nu1(1)+4*nu1(2)-nu1(3))/(2*s_elem1(cl2mem11));
dnuds1(2:Nm1-1) = (nu1(3:Nm1)-nu1(1:Nm1-2))./(s_elem1(cl2mem11:cl2mem12-2)+s_elem1(cl2mem11+1:cl2mem12-1));
dnuds1(Nm1) = (3*nu1(Nm1)-4*nu1(Nm1-1)+nu1(Nm1-2))/(2*s_elem1(cl2mem12-1));

% compute d2kappads
d2kappads1 = zeros(Nm1,1);
d2kappads1(1) = ...
    (2*kappa1(cl2mem11)-5*kappa1(cl2mem11+1)+4*kappa1(cl2mem11+2)-kappa1(cl2mem11+3))...
    ./(s_elem1(cl2mem11)^2);
d2kappads1(2:Nm1-1) = ...
    (kappa1(cl2mem11:cl2mem12-2)+kappa1(cl2mem11+2:cl2mem12)-2*kappa1(cl2mem11+1:cl2mem12-1))...
    ./(((s_elem1(cl2mem11:cl2mem12-2)+s_elem1(cl2mem11+1:cl2mem12-1))/2).^2);
d2kappads1(Nm1) = ...
    (2*kappa1(cl2mem12)-5*kappa1(cl2mem12-1)+4*kappa1(cl2mem12-2)-kappa1(cl2mem12-3))...
    ./(s_elem1(cl2mem12-1)^2);


d2kappads1(Nm1-1) = 2*d2kappads1(Nm1-2)-d2kappads1(Nm1-1);
d2kappads1(Nm1) = 2*d2kappads1(Nm1-1)-d2kappads1(Nm1);

% p1 = plot(x1(cl2mem11:cl2mem12),-cb1/Ca1*d2kappads1,'-','LineWidth',2.0,'color',c1);
% hold on;
% p2 = plot(x1(cl2mem11:cl2mem12),-cb1/2/Ca1*kappa1(cl2mem11:cl2mem12).^3,'--','LineWidth',2.0,'color',c2);
% p3 = plot(x1(cl2mem11:cl2mem12),1/Ca1*(gamma11-nu1(1:Nm1)).*kappa1(cl2mem11:cl2mem12),':','LineWidth',2.0,'color',c3);
% % p4 = plot(x1(cl2mem11:cl2mem12),1/Ca1*dnuds1,'-.','LineWidth',2.0,'color',c4);
% p5 = plot(x1(cl2mem11:cl2mem12),pressure1(mem2node1),'--k','LineWidth',2.0);
% xlabel('$x$','Interpreter','latex');
% hh = legend([p1,p2,p3,p5],{'$-\frac{c_b}{Ca}\Delta_s\kappa$','$-\frac{c_b}{2Ca}\kappa^3$',...
%     '$\frac{1}{Ca}(\gamma_1-\nu)\kappa$','$p|_{\Sigma_1}$'});
% set(hh,'Interpreter','latex','FontName','Times New Roman','FontSize',32,'FontWeight','normal','Location','southeast');
% set(gca,'FontSize',32);
% xlim([x1(cl2mem11-5), x1(cl2mem12+5)]);
% hold off;

ck = [0 0 0]/256;
fig = figure;
set(fig,'defaultAxesColorOrder',[ck; c2]);
yyaxis left
p1 = plot(x1(cl2mem11:cl2mem12),cb1/Ca1*d2kappads1+cb1/2/Ca1*kappa1(cl2mem11:cl2mem12).^3,'-','LineWidth',2.0,'color',ck);
hold on;
% p5 = plot(x1(cl2mem11:cl2mem12),pressure1(mem2node1),'--','LineWidth',2.0,'color',ck);
% p4 = plot(x1(cl2mem11:cl2mem12),1/Ca1*dnuds1,'-.','LineWidth',2.0,'color',c2);

yyaxis right
% p2 = plot(x1(cl2mem11:cl2mem12),-cb1/2/Ca1*kappa1(cl2mem11:cl2mem12).^3,'-','LineWidth',2.0,'color',c2);
p3 = plot(x1(cl2mem11:cl2mem12),1/Ca1*(gamma11-nu1(1:Nm1)).*kappa1(cl2mem11:cl2mem12),'-','LineWidth',2.0,'color',c2);
% p4 = plot(x1(cl2mem11:cl2mem12),1/Ca1*dnuds1,'-.','LineWidth',2.0,'color',c2);
xlabel('$x$','Interpreter','latex');
% hh = legend([p1,p5,p2,p3,p4],{'$-\frac{c_b}{Ca}\Delta_s\kappa$','$p|_{\Sigma_1}$','$-\frac{c_b}{2Ca}\kappa^3$','$\frac{1}{Ca}(\gamma_1-\nu)\kappa$', '\unboldmath{$\frac{1}{Ca}$}\boldmath{$\tau\cdot\nabla$}\unboldmath{$_s\nu$}'});
% hh = legend([p1,p5,p2,p3],{'$-\frac{c_b}{Ca}\Delta_s\kappa$','$p|_{\Sigma_1}$',...
%             '$-\frac{c_b}{2Ca}\kappa^3$','$\frac{1}{Ca}(\gamma_1-\nu)\kappa$'});
hh = legend([p1,p3],{'$\frac{c_b}{Ca}(\Delta_s\kappa+\frac{1}{2}\kappa^3$','$\frac{1}{Ca}(\gamma_1-\nu)\kappa$'});
        
set(hh,'Interpreter','latex','FontName','Times New Roman','FontSize',32,'FontWeight','normal','Location','southeast');
set(gca,'FontSize',32);
xlim([x1(cl2mem11-5), x1(cl2mem12+5)]);
hold off;
