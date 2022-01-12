clear;
clc;
addpath(genpath('.'));
dbstop if error

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end

load('record_wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1_ind_17.mat');
nu1 = sln.nu;
pmem1 = sln.pmemold; Nmem1 = para.Nmem; Nm1 = para.Nm; x1 = pmem1(:,1);
cl2mem11 = geom.cl2mem(1); cl2mem12 = geom.cl2mem(2);


load('record_dewetting_20_10_gap_4_xl_2_gamma2_0_3_cb_5e3_Ca_0_1_ind_14.mat');
nu2 = sln.nu;
pmem2 = sln.pmemold; Nmem2 = para.Nmem; Nm2 = para.Nm; x2 = pmem2(:,1);
cl2mem21 = geom.cl2mem(1); cl2mem22 = geom.cl2mem(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot nu vs x

p2 = plot(x2(cl2mem21:cl2mem22),nu2(1:Nm2),'-','LineWidth',2.0,'color',c2);
hold on;
plot([x2(1) x2(cl2mem21)], [nu2(Nm2+1) nu2(Nm2+1)],'-','LineWidth',2.0,'color',c2);
plot([x2(cl2mem22) x2(Nmem2)], [0 0],'-','LineWidth',2.0,'color',c2);

p1 = plot(x1(cl2mem11:cl2mem12),nu1(1:Nm1),'--','LineWidth',2.0,'color',c1);
plot([x1(1) x1(cl2mem11)], [nu1(Nm1+1) nu1(Nm1+1)],'--','LineWidth',2.0,'color',c1);
plot([x1(cl2mem12) x1(Nmem1)], [0 0],'--','LineWidth',2.0,'color',c1);

hh = legend([p1,p2],{'wetting: $\cos\theta_Y=0.5$','dewetting: $\cos\theta_Y=-0.7$'});
set(hh,'Interpreter','latex','FontName','Times New Roman','FontSize',24,'FontWeight','normal','Location','southeast');

xlabel('$x$','Interpreter','latex');
ylabel('$\nu$','Interpreter','latex');
set(gca,'FontSize',32);
hold off;