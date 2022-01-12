clear;
clc;
addpath(genpath('.'));
dbstop if error

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end

load('record_wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1_ind_1.mat'); 
cb1 = para.cb; kappa1 = sln.kmem;
pmem1 = sln.pmemold; Nmem1 = para.Nmem; x1 = pmem1(:,1); y1 = pmem1(:,2);
cl2mem11 = geom.cl2mem(1); cl2mem12 = geom.cl2mem(2);
s_elem1 = sqrt((x1(2:end)-x1(1:end-1)).^2+(y1(2:end)-y1(1:end-1)).^2);

dkappads1 = zeros(Nmem1+2,1);
dkappads1(1) = (-3*kappa1(1)+4*kappa1(2)-kappa1(3))/(2*s_elem1(1));
dkappads1(2:cl2mem11-1) = (kappa1(3:cl2mem11)-kappa1(1:cl2mem11-2))./(s_elem1(1:cl2mem11-2)+s_elem1(2:cl2mem11-1));
dkappads1(cl2mem11) = (3*kappa1(cl2mem11)-4*kappa1(cl2mem11-1)+kappa1(cl2mem11-2))/(2*s_elem1(cl2mem11-1));

dkappads1(Nmem1+1) = (-3*kappa1(cl2mem11)+4*kappa1(cl2mem11+1)-kappa1(cl2mem11+2))/(2*s_elem1(cl2mem11));
dkappads1(cl2mem11+1:cl2mem12-1) = (kappa1(cl2mem11+2:cl2mem12)-kappa1(cl2mem11:cl2mem12-2))./(s_elem1(cl2mem11:cl2mem12-2)+s_elem1(cl2mem11+1:cl2mem12-1));
dkappads1(Nmem1+2) = (3*kappa1(cl2mem12)-4*kappa1(cl2mem12-1)+kappa1(cl2mem12-2))/(2*s_elem1(cl2mem12-1));

dkappads1(cl2mem12) = (-3*kappa1(cl2mem12)+4*kappa1(cl2mem12+1)-kappa1(cl2mem12+2))/(2*s_elem1(cl2mem12));
dkappads1(cl2mem12+1:Nmem1-1) = (kappa1(cl2mem12+2:Nmem1)-kappa1(cl2mem12:Nmem1-2))./(s_elem1(cl2mem12:Nmem1-2)+s_elem1(cl2mem12+1:Nmem1-1));
dkappads1(Nmem1) = (3*kappa1(Nmem1)-4*kappa1(Nmem1-1)+kappa1(Nmem1-2))/(2*s_elem1(Nmem1-1));

% ------------------------------------------------------------------------------------------------------
load('record_dewetting_20_10_gap_4_xl_2_gamma2_0_3_cb_5e3_Ca_0_1_ind_1.mat');

cb2 = para.cb; kappa2 = sln.kmem;
pmem2 = sln.pmemold; Nmem2 = para.Nmem; x2 = pmem2(:,1); y2 = pmem2(:,2);
cl2mem21 = geom.cl2mem(1); cl2mem22 = geom.cl2mem(2);
s_elem2 = sqrt((x2(2:end)-x2(1:end-1)).^2+(y2(2:end)-y2(1:end-1)).^2);

dkappads2 = zeros(Nmem2+2,1);
dkappads2(1) = (-3*kappa2(1)+4*kappa2(2)-kappa2(3))/(2*s_elem2(1));
dkappads2(2:cl2mem21-1) = (kappa2(3:cl2mem21)-kappa2(1:cl2mem21-2))./(s_elem2(1:cl2mem21-2)+s_elem2(2:cl2mem21-1));
dkappads2(cl2mem21) = (3*kappa2(cl2mem21)-4*kappa2(cl2mem21-1)+kappa2(cl2mem21-2))/(2*s_elem2(cl2mem21-1));

dkappads2(Nmem2+1) = (-3*kappa2(cl2mem21)+4*kappa2(cl2mem21+1)-kappa2(cl2mem21+2))/(2*s_elem2(cl2mem21));
dkappads2(cl2mem21+1:cl2mem22-1) = (kappa2(cl2mem21+2:cl2mem22)-kappa2(cl2mem21:cl2mem22-2))./(s_elem2(cl2mem21:cl2mem22-2)+s_elem2(cl2mem21+1:cl2mem22-1));
dkappads2(Nmem2+2) = (3*kappa2(cl2mem22)-4*kappa2(cl2mem22-1)+kappa2(cl2mem22-2))/(2*s_elem2(cl2mem22-1));

dkappads2(cl2mem22) = (-3*kappa2(cl2mem22)+4*kappa2(cl2mem22+1)-kappa2(cl2mem22+2))/(2*s_elem2(cl2mem22));
dkappads2(cl2mem22+1:Nmem2-1) = (kappa2(cl2mem22+2:Nmem2)-kappa2(cl2mem22:Nmem2-2))./(s_elem2(cl2mem22:Nmem2-2)+s_elem2(cl2mem22+1:Nmem2-1));
dkappads2(Nmem2) = (3*kappa2(Nmem2)-4*kappa2(Nmem2-1)+kappa2(Nmem2-2))/(2*s_elem2(Nmem2-1));
% ----------------------------------------------------------------------------------------------------


p3 = plot(x2(1:Nmem2),kappa2(1:Nmem2),'-','LineWidth',2.0,'color',c2);
hold on;
plot([x2(cl2mem21) x2(cl2mem22)],[kappa2(cl2mem21) kappa2(cl2mem22)],'sk','LineWidth',2.0);
p4 = plot(x1(1:Nmem1),kappa1(1:Nmem1),'-.','LineWidth',2.0,'color',c1);plot(x2(1:Nmem2),kappa2(1:Nmem2),'-','LineWidth',2.0,'color',c2);
plot([x1(cl2mem11) x1(cl2mem12)],[kappa1(cl2mem11) kappa1(cl2mem12)],'ok','LineWidth',2.0);
hh = legend([p4,p3],{'wetting: $\cos\theta_Y=0.5$','dewetting: $\cos\theta_Y=-0.7$'});
set(hh,'Interpreter','latex','FontName','Times New Roman','FontSize',32,'FontWeight','normal','Location','southeast');

xlabel('$x$','Interpreter','latex');
ylabel('$\kappa$','Interpreter','latex');
set(gca,'FontSize',32);
hold off;

% p1 = plot(x1(1:cl2mem11),dkappads1(1:cl2mem11),'-.','LineWidth',2.0,'color',c1);
% hold on;
% plot(x1(cl2mem11:cl2mem12),[dkappads1(Nmem1+1);dkappads1(cl2mem11+1:cl2mem12-1);dkappads1(Nmem1+2)],'-.','LineWidth',2.0,'color',c1);
% plot(x1(cl2mem12:Nmem1),dkappads1(cl2mem12:Nmem1),'-.','LineWidth',2.0,'color',c1);
% 
% p2 = plot(x2(1:cl2mem21),dkappads2(1:cl2mem21),'-','LineWidth',2.0,'color',c2);
% plot(x2(cl2mem21:cl2mem22),[dkappads2(Nmem2+1);dkappads2(cl2mem21+1:cl2mem22-1);dkappads2(Nmem2+2)],'-','LineWidth',2.0,'color',c2);
% plot(x2(cl2mem22:Nmem2),dkappads2(cl2mem22:Nmem2),'-','LineWidth',2.0,'color',c2);
% 
% hh = legend([p1,p2],{'wetting: $\cos\theta_Y=0.5$','dewetting: $\cos\theta_Y=-0.7$'});
% set(hh,'Interpreter','latex','FontName','Times New Roman','FontSize',32,'FontWeight','normal','Location','southeast');
% 
% xlabel('$x$','Interpreter','latex');
% ylabel('\unboldmath{$c_b$}\boldmath{$\tau\cdot\nabla$}\unboldmath{$_s\kappa$}','Interpreter','latex');
% set(gca,'FontSize',32);
% hold off;


p1 = plot(x1(1:cl2mem11),cb1*dkappads1(1:cl2mem11),'-.','LineWidth',2.0,'color',c1);
hold on;
plot(x1(cl2mem11:cl2mem12),cb1*[dkappads1(Nmem1+1);dkappads1(cl2mem11+1:cl2mem12-1);dkappads1(Nmem1+2)],'-.','LineWidth',2.0,'color',c1);
plot(x1(cl2mem12:Nmem1),cb1*dkappads1(cl2mem12:Nmem1),'-.','LineWidth',2.0,'color',c1);

p2 = plot(x2(1:cl2mem21),cb2*dkappads2(1:cl2mem21),'-','LineWidth',2.0,'color',c2);
plot(x2(cl2mem21:cl2mem22),cb2*[dkappads2(Nmem2+1);dkappads2(cl2mem21+1:cl2mem22-1);dkappads2(Nmem2+2)],'-','LineWidth',2.0,'color',c2);
plot(x2(cl2mem22:Nmem2),cb2*dkappads2(cl2mem22:Nmem2),'-','LineWidth',2.0,'color',c2);

hh = legend([p1,p2],{'wetting: $\cos\theta_Y=0.5,\ \sin\theta_Y=0.866$','dewetting: $\cos\theta_Y=-0.7,\ \sin\theta_Y=0.714$'});
set(hh,'Interpreter','latex','FontName','Times New Roman','FontSize',32,'FontWeight','normal','Location','southeast');

xlabel('$x$','Interpreter','latex');
ylabel('\unboldmath{$c_b$}\boldmath{$\tau\cdot\nabla$}\unboldmath{$_s\kappa$}','Interpreter','latex');
set(gca,'FontSize',32);
hold off;

cb1*(dkappads1(Nmem1+1)-dkappads1(cl2mem11))
cb1*(dkappads1(Nmem1+2)-dkappads1(cl2mem12))

cb2*(dkappads2(Nmem2+1)-dkappads2(cl2mem21))
cb2*(dkappads2(Nmem2+2)-dkappads2(cl2mem22))