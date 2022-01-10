% close all;
clear;
% clc;
addpath(genpath('.'));
dbstop if error


for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end 

% load('rectangle_wetting_data128_cb1e_1_gamma1_0_5_reprod_Ca_0_2.mat');
% load('rectangle_dewetting_data256_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');
load('record_curvature_rect_dewetting_data128.mat');
% load('record_curvature_rect_wetting_data128.mat');

ind = floor(0.05/para.dt)+1; Nffi = para.Nffi; Nmem = hstry.Nmems(ind);
pmem = reshape(hstry.pmem_track(ind-1,1:2*Nmem-1,1:2),2*Nmem-1,2);
kmem = reshape(hstry.kmem_track(ind,1:2*Nmem-1),2*Nmem-1,1);
    pmem2 = zeros(2*Nmem-1,2);
    pmem2(1:2:2*Nmem-1,1:2) = pmem(1:Nmem,1:2);
    pmem2(2:2:2*Nmem-2,1:2) = pmem(Nmem+1:2*Nmem-1,1:2);
    kmem2 = zeros(2*Nmem-1,1);
    kmem2(1:2:2*Nmem-1) = kmem(1:Nmem);
    kmem2(2:2:2*Nmem-2) = kmem(Nmem+1:2*Nmem-1);
    
plot(pmem2(:,1),kmem2,'-','LineWidth',2.0,'color',c1);
hold on;

ind = floor(0.2/para.dt)+1;
pmem = reshape(hstry.pmem_track(ind-1,1:2*Nmem-1,1:2),2*Nmem-1,2);
kmem = reshape(hstry.kmem_track(ind,1:2*Nmem-1),2*Nmem-1,1);
    pmem2 = zeros(2*Nmem-1,2);
    pmem2(1:2:2*Nmem-1,1:2) = pmem(1:Nmem,1:2);
    pmem2(2:2:2*Nmem-2,1:2) = pmem(Nmem+1:2*Nmem-1,1:2);
    kmem2 = zeros(2*Nmem-1,1);
    kmem2(1:2:2*Nmem-1) = kmem(1:Nmem);
    kmem2(2:2:2*Nmem-2) = kmem(Nmem+1:2*Nmem-1);
    
plot(pmem2(:,1),kmem2,'--','LineWidth',2.0,'color',c2);

ind = floor(0.5/para.dt)+1;
pmem = reshape(hstry.pmem_track(ind-1,1:2*Nmem-1,1:2),2*Nmem-1,2);
kmem = reshape(hstry.kmem_track(ind,1:2*Nmem-1),2*Nmem-1,1);
    pmem2 = zeros(2*Nmem-1,2);
    pmem2(1:2:2*Nmem-1,1:2) = pmem(1:Nmem,1:2);
    pmem2(2:2:2*Nmem-2,1:2) = pmem(Nmem+1:2*Nmem-1,1:2);
    kmem2 = zeros(2*Nmem-1,1);
    kmem2(1:2:2*Nmem-1) = kmem(1:Nmem);
    kmem2(2:2:2*Nmem-2) = kmem(Nmem+1:2*Nmem-1);
    
plot(pmem2(:,1),kmem2,'-.','LineWidth',2.0,'color',c3);

ind = floor(1/para.dt)+1;
pmem = reshape(hstry.pmem_track(ind-1,1:2*Nmem-1,1:2),2*Nmem-1,2);
kmem = reshape(hstry.kmem_track(ind,1:2*Nmem-1),2*Nmem-1,1);
    pmem2 = zeros(2*Nmem-1,2);
    pmem2(1:2:2*Nmem-1,1:2) = pmem(1:Nmem,1:2);
    pmem2(2:2:2*Nmem-2,1:2) = pmem(Nmem+1:2*Nmem-1,1:2);
    kmem2 = zeros(2*Nmem-1,1);
    kmem2(1:2:2*Nmem-1) = kmem(1:Nmem);
    kmem2(2:2:2*Nmem-2) = kmem(Nmem+1:2*Nmem-1);
    
plot(pmem2(:,1),kmem2,'-','LineWidth',2.0,'color',c4);

ind = floor(2/para.dt)+1;
pmem = reshape(hstry.pmem_track(ind-1,1:2*Nmem-1,1:2),2*Nmem-1,2);
kmem = reshape(hstry.kmem_track(ind,1:2*Nmem-1),2*Nmem-1,1);
    pmem2 = zeros(2*Nmem-1,2);
    pmem2(1:2:2*Nmem-1,1:2) = pmem(1:Nmem,1:2);
    pmem2(2:2:2*Nmem-2,1:2) = pmem(Nmem+1:2*Nmem-1,1:2);
    kmem2 = zeros(2*Nmem-1,1);
    kmem2(1:2:2*Nmem-1) = kmem(1:Nmem);
    kmem2(2:2:2*Nmem-2) = kmem(Nmem+1:2*Nmem-1);
    
plot(pmem2(:,1),kmem2,':','LineWidth',2.0,'color',c5);

leg = legend('$t=0.05$','$t=0.2$','$t=0.5$','$t=1$','$t=2$');
set(leg,'Interpreter','latex'); 
xlabel('$x$','Interpreter','latex');
ylabel('$\kappa$','Interpreter','latex');
set(gca,'FontSize',32);
axis equal;
% ylim([-0.5,1.1]);
ylim([-0.85,1.1]);
xlim([-1,1]);
hold off;
