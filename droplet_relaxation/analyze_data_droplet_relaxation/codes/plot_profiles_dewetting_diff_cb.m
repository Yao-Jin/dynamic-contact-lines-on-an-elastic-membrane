% close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time = 2; lw = 2;

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end 

load('rectangle_dewetting_data128_cb1e2_gamma2_0_5_noprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
p1 = plot(hstry.pmem_track(ind,1:hstry.Nmems(ind),1),...
     hstry.pmem_track(ind,1:hstry.Nmems(ind),2),'-k','LineWidth',lw);
hold on;
plot(hstry.pffi_track(ind,1:para.Nffi,1),...
     hstry.pffi_track(ind,1:para.Nffi,2),'-k','LineWidth',lw);

load('rectangle_dewetting_data256_cb1e_3_gamma2_0_5_reprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
p2 = plot(hstry.pmem_track(ind,1:hstry.Nmems(ind),1),...
     hstry.pmem_track(ind,1:hstry.Nmems(ind),2),'-','LineWidth',lw,'color',c3);
plot(hstry.pffi_track(ind,1:para.Nffi,1),...
     hstry.pffi_track(ind,1:para.Nffi,2),'-','LineWidth',lw,'color',c3);
 
% load('rectangle_dewetting_data128_cb1e0_gamma2_0_5_noprod_Ca_0_2.mat');
% ind = floor(Time/para.dt)+1;
% p3 = plot(hstry.pmem_track(ind,1:hstry.Nmems(ind),1),...
%      hstry.pmem_track(ind,1:hstry.Nmems(ind),2),'--','LineWidth',lw,'color',c3);
% plot(hstry.pffi_track(ind,1:para.Nffi,1),...
%      hstry.pffi_track(ind,1:para.Nffi,2),'--','LineWidth',lw,'color',c3);
 
load('rectangle_dewetting_data128_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
p4 = plot(hstry.pmem_track(ind,1:hstry.Nmems(ind),1),...
     hstry.pmem_track(ind,1:hstry.Nmems(ind),2),'--','LineWidth',lw,'color',c2);
plot(hstry.pffi_track(ind,1:para.Nffi,1),...
     hstry.pffi_track(ind,1:para.Nffi,2),'--','LineWidth',lw,'color',c2);
 
load('rectangle_dewetting_data128_cb1e_2_gamma2_0_5_noprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
p5 = plot(hstry.pmem_track(ind,1:hstry.Nmems(ind),1),...
     hstry.pmem_track(ind,1:hstry.Nmems(ind),2),'-.','LineWidth',lw,'color',c1);
plot(hstry.pffi_track(ind,1:para.Nffi,1),...
     hstry.pffi_track(ind,1:para.Nffi,2),'-.','LineWidth',lw,'color',c1);
 
leg = legend([p1,p4,p5,p2],{'$c_b=100$','$c_b=0.1$','$c_b=0.01$','$c_b=0.001$'});
set(leg,'Interpreter','latex'); 
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
set(gca,'FontSize',32);
axis equal;
ylim([-0.3,0.7]); 
xlim([-1,1]);
hold off;
