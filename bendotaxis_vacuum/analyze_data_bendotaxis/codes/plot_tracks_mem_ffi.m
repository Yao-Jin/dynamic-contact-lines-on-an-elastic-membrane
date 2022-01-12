% close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time = 3; lw = 3;
% disp(['Time =',num2str(Time)]);
% disp('cb:1e_1 wetting:');

load('testconv_10_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');
ind = floor(Time/para.dt)+1; Bx = para.Bx; By = para.By;
p1 = plot(hstry.pmem_track(ind,1:hstry.Nmems(ind),1),hstry.pmem_track(ind,1:hstry.Nmems(ind),2)-By,':k','LineWidth',lw);
hold on;
plot(hstry.pffil_track(ind,1:hstry.Nffils(ind),1),hstry.pffil_track(ind,1:hstry.Nffils(ind),2)-By,':k','LineWidth',lw);
plot(hstry.pffir_track(ind,1:hstry.Nffirs(ind),1),hstry.pffir_track(ind,1:hstry.Nffirs(ind),2)-By,':k','LineWidth',lw);

load('testconv_20_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');
ind = floor(Time/para.dt)+1;
p2 = plot(hstry.pmem_track(ind,1:hstry.Nmems(ind),1),hstry.pmem_track(ind,1:hstry.Nmems(ind),2)-By,'--','LineWidth',lw,'color',c2);
plot(hstry.pffil_track(ind,1:hstry.Nffils(ind),1),hstry.pffil_track(ind,1:hstry.Nffils(ind),2)-By,'--','LineWidth',lw,'color',c2);
plot(hstry.pffir_track(ind,1:hstry.Nffirs(ind),1),hstry.pffir_track(ind,1:hstry.Nffirs(ind),2)-By,'--','LineWidth',lw,'color',c2);

load('testconv_40_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');
ind = floor(Time/para.dt)+1;
p3 = plot(hstry.pmem_track(ind,1:hstry.Nmems(ind),1),hstry.pmem_track(ind,1:hstry.Nmems(ind),2)-By,'-.','LineWidth',lw,'color',c3);
plot(hstry.pffil_track(ind,1:hstry.Nffils(ind),1),hstry.pffil_track(ind,1:hstry.Nffils(ind),2)-By,'-.','LineWidth',lw,'color',c3);
plot(hstry.pffir_track(ind,1:hstry.Nffirs(ind),1),hstry.pffir_track(ind,1:hstry.Nffirs(ind),2)-By,'-.','LineWidth',lw,'color',c3);

load('testconv_80_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');
ind = floor(Time/para.dt)+1;
p4 = plot(hstry.pmem_track(ind,1:hstry.Nmems(ind),1),hstry.pmem_track(ind,1:hstry.Nmems(ind),2)-By,'-','LineWidth',lw,'color',c1);
plot(hstry.pffil_track(ind,1:hstry.Nffils(ind),1),hstry.pffil_track(ind,1:hstry.Nffils(ind),2)-By,'-','LineWidth',lw,'color',c1);
plot(hstry.pffir_track(ind,1:hstry.Nffirs(ind),1),hstry.pffir_track(ind,1:hstry.Nffirs(ind),2)-By,'-','LineWidth',lw,'color',c1);

plot([0,Bx],[0 0],'-k','LineWidth',1.0);

leg = legend([p1,p2,p3,p4],{'$\mathcal{R}=10$', '$\mathcal{R}=20$', '$\mathcal{R}=40$', '$\mathcal{R}=80$'});
set(leg,'Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
set(gca,'FontSize',24);
axis equal;
% ylim([-0.5,1]); 
% xlim([-1,1]);
hold off;
