clear;
clc;
addpath(genpath('.'));
dbstop if error

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end

% load('record_wetting_gamma1_0_5_cb_5e3_t_200.mat');
% load('record_dewetting_gamma2_0_3_cb_5e3_t_600.mat');

%1 4 9 17
% load('record_wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1_ind_1.mat'); 
%1 4 9 14
load('record_dewetting_20_10_gap_4_xl_2_gamma2_0_3_cb_5e3_Ca_0_1_ind_1.mat');
disp(['Time = ',num2str((hstry.itr-1)*para.dt)]);

cb = para.cb; gamma1 = para.gamma1; Ca = para.Ca;
nu1 = sln.nu; kappa1 = sln.kmem;
pmem1 = sln.pmemold; Nmem1 = para.Nmem; Nm1 = para.Nm; x1 = pmem1(1:Nmem1,1); y1 = pmem1(1:Nmem1,2);
cl2mem11 = geom.cl2mem(1); cl2mem12 = geom.cl2mem(2);
s_elem1 = sqrt((x1(2:end)-x1(1:end-1)).^2+(y1(2:end)-y1(1:end-1)).^2);
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
% d2kappads1(Nm1-1) = 2*d2kappads1(Nm1-2)-d2kappads1(Nm1-1);
% d2kappads1(Nm1) = 2*d2kappads1(Nm1-1)-d2kappads1(Nm1);
%------------------------------------------------------------------------------------------------

%1 4 9 17
% load('record_wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1_ind_4.mat'); 
%1 4 9 14
load('record_dewetting_20_10_gap_4_xl_2_gamma2_0_3_cb_5e3_Ca_0_1_ind_4.mat');
disp(['Time = ',num2str((hstry.itr-1)*para.dt)]);

nu2 = sln.nu; kappa2 = sln.kmem;
pmem2 = sln.pmemold; Nmem2 = para.Nmem; Nm2 = para.Nm; x2 = pmem2(1:Nmem2,1); y2 = pmem2(1:Nmem2,2);
cl2mem21 = geom.cl2mem(1); cl2mem22 = geom.cl2mem(2);
s_elem2 = sqrt((x2(2:end)-x2(1:end-1)).^2+(y2(2:end)-y2(1:end-1)).^2);
% compute d2kappads
d2kappads2 = zeros(Nm2,1);
d2kappads2(1) = ...
    (2*kappa2(cl2mem21)-5*kappa2(cl2mem21+1)+4*kappa2(cl2mem21+2)-kappa2(cl2mem21+3))...
    ./(s_elem2(cl2mem21)^2);
d2kappads2(2:Nm2-1) = ...
    (kappa2(cl2mem21:cl2mem22-2)+kappa2(cl2mem21+2:cl2mem22)-2*kappa2(cl2mem21+1:cl2mem22-1))...
    ./(((s_elem2(cl2mem21:cl2mem22-2)+s_elem2(cl2mem21+1:cl2mem22-1))/2).^2);
d2kappads2(Nm2) = ...
    (2*kappa2(cl2mem22)-5*kappa2(cl2mem22-1)+4*kappa2(cl2mem22-2)-kappa2(cl2mem22-3))...
    ./(s_elem2(cl2mem22-1)^2);
%------------------------------------------------------------------------------------------------

%1 4 9 17
% load('record_wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1_ind_9.mat'); 
%1 4 9 14
load('record_dewetting_20_10_gap_4_xl_2_gamma2_0_3_cb_5e3_Ca_0_1_ind_9.mat');
disp(['Time = ',num2str((hstry.itr-1)*para.dt)]);

nu3 = sln.nu; kappa3 = sln.kmem;
pmem3 = sln.pmemold; Nmem3 = para.Nmem; Nm3 = para.Nm; x3 = pmem3(1:Nmem3,1); y3 = pmem3(1:Nmem3,2);
cl2mem31 = geom.cl2mem(1); cl2mem32 = geom.cl2mem(2);
s_elem3 = sqrt((x3(2:end)-x3(1:end-1)).^2+(y3(2:end)-y3(1:end-1)).^2);
% compute d2kappads
d2kappads3 = zeros(Nm3,1);
d2kappads3(1) = ...
    (2*kappa3(cl2mem31)-5*kappa3(cl2mem31+1)+4*kappa3(cl2mem31+2)-kappa3(cl2mem31+3))...
    ./(s_elem3(cl2mem31)^2);
d2kappads3(2:Nm3-1) = ...
    (kappa3(cl2mem31:cl2mem32-2)+kappa3(cl2mem31+2:cl2mem32)-2*kappa3(cl2mem31+1:cl2mem32-1))...
    ./(((s_elem3(cl2mem31:cl2mem32-2)+s_elem3(cl2mem31+1:cl2mem32-1))/2).^2);
d2kappads3(Nm3) = ...
    (2*kappa3(cl2mem32)-5*kappa3(cl2mem32-1)+4*kappa3(cl2mem32-2)-kappa3(cl2mem32-3))...
    ./(s_elem3(cl2mem32-1)^2);
%------------------------------------------------------------------------------------------------

%1 4 9 17
% load('record_wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1_ind_17.mat'); 
%1 4 9 14
load('record_dewetting_20_10_gap_4_xl_2_gamma2_0_3_cb_5e3_Ca_0_1_ind_14.mat');
disp(['Time = ',num2str((hstry.itr-1)*para.dt)]);

nu4 = sln.nu; kappa4 = sln.kmem;
pmem4 = sln.pmemold; Nmem4 = para.Nmem; Nm4 = para.Nm; x4 = pmem4(1:Nmem4,1); y4 = pmem4(1:Nmem4,2);
cl2mem41 = geom.cl2mem(1); cl2mem42 = geom.cl2mem(2);
s_elem4 = sqrt((x4(2:end)-x4(1:end-1)).^2+(y4(2:end)-y4(1:end-1)).^2);
% compute d2kappads
d2kappads4 = zeros(Nm4,1);
d2kappads4(1) = ...
    (2*kappa4(cl2mem41)-5*kappa4(cl2mem41+1)+4*kappa4(cl2mem41+2)-kappa4(cl2mem41+3))...
    ./(s_elem4(cl2mem41)^2);
d2kappads4(2:Nm4-1) = ...
    (kappa4(cl2mem41:cl2mem42-2)+kappa4(cl2mem41+2:cl2mem42)-2*kappa4(cl2mem41+1:cl2mem42-1))...
    ./(((s_elem4(cl2mem41:cl2mem42-2)+s_elem4(cl2mem41+1:cl2mem42-1))/2).^2);
d2kappads4(Nm4) = ...
    (2*kappa4(cl2mem42)-5*kappa4(cl2mem42-1)+4*kappa4(cl2mem42-2)-kappa4(cl2mem42-3))...
    ./(s_elem4(cl2mem42-1)^2);


% yyaxis left
p1 = plot(x1(cl2mem11:cl2mem12),cb/Ca*d2kappads1+cb/2/Ca*kappa1(cl2mem11:cl2mem12).^3,':','LineWidth',2.0,'color',c1);
hold on;
p2 = plot(x2(cl2mem21:cl2mem22),cb/Ca*d2kappads2+cb/2/Ca*kappa2(cl2mem21:cl2mem22).^3,'--','LineWidth',2.0,'color',c2);
p3 = plot(x3(cl2mem31:cl2mem32),cb/Ca*d2kappads3+cb/2/Ca*kappa3(cl2mem31:cl2mem32).^3,'-.','LineWidth',2.0,'color',c3);
p4 = plot(x4(cl2mem41:cl2mem42),cb/Ca*d2kappads4+cb/2/Ca*kappa4(cl2mem41:cl2mem42).^3,'-','LineWidth',2.0,'color',c4);
xlabel('$x$','Interpreter','latex');
ylabel('$\frac{c_b}{Ca}(\Delta_s\kappa+\frac{1}{2}\kappa^3)$','Interpreter','latex');
hh = legend([p1,p2,p3,p4],{'$t=0.05$','$t=0.3$','$t=1.0$','$t=15.0$'});
set(hh,'Interpreter','latex','FontName','Times New Roman','FontSize',32,'FontWeight','normal','Location','northwest');
set(gca,'FontSize',32);
xlim([x1(cl2mem11-5), x4(cl2mem42+5)]);
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p5 = plot(x1(cl2mem11:cl2mem12),1/Ca*(nu1(1:Nm1)-gamma1).*kappa1(cl2mem11:cl2mem12),':','LineWidth',2.0,'color',c1);
hold on;
p6 = plot(x2(cl2mem21:cl2mem22),1/Ca*(nu2(1:Nm2)-gamma1).*kappa2(cl2mem21:cl2mem22),'--','LineWidth',2.0,'color',c2);
p7 = plot(x3(cl2mem31:cl2mem32),1/Ca*(nu3(1:Nm3)-gamma1).*kappa3(cl2mem31:cl2mem32),'-.','LineWidth',2.0,'color',c3);
p8 = plot(x4(cl2mem41:cl2mem42),1/Ca*(nu4(1:Nm4)-gamma1).*kappa4(cl2mem41:cl2mem42),'-','LineWidth',2.0,'color',c4);
xlabel('$x$','Interpreter','latex');
ylabel('$\frac{1}{Ca}(\nu-\gamma_1)\kappa$','Interpreter','latex');
hh = legend([p5,p6,p7,p8],{'$t=0.05$','$t=0.3$','$t=1.0$','$t=15.0$'});
set(hh,'Interpreter','latex','FontName','Times New Roman','FontSize',32,'FontWeight','normal','Location','northeast');
set(gca,'FontSize',32);
xlim([x1(cl2mem11-5), x4(cl2mem42+5)]);
hold off;