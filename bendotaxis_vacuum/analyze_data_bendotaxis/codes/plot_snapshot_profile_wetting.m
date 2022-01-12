% close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error


% cos(thetaY) = (gamma2-gamma1)/gamma3; gamma3 = 1;
% wetting: gamma2 = 1;
% dewetting: gamma1 = 1;

for readdata = 1:1
load('wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1.mat');
lim1 = hstry.itr; pmem1 = hstry.pmem_track; pffil1 = hstry.pffil_track; pffir1 = hstry.pffir_track;
Nmem1 = hstry.Nmems; Nffil1 = hstry.Nffils; Nffir1 = hstry.Nffirs;
pmem1(:,:,2) = pmem1(:,:,2)-para.By; pffil1(:,:,2) = pffil1(:,:,2)-para.By; pffir1(:,:,2) = pffir1(:,:,2)-para.By;

load('wetting_20_10_gap_4_xl_2_gamma1_0_5_cb_5e3_Ca_0_1.mat');
lim2 = hstry.itr; pmem2 = hstry.pmem_track; pffil2 = hstry.pffil_track; pffir2 = hstry.pffir_track;
Nmem2 = hstry.Nmems; Nffil2 = hstry.Nffils; Nffir2 = hstry.Nffirs;
pmem2(:,:,2) = pmem2(:,:,2)-para.By; pffil2(:,:,2) = pffil2(:,:,2)-para.By; pffir2(:,:,2) = pffir2(:,:,2)-para.By;

load('wetting_20_10_gap_4_xl_2_gamma1_0_7_cb_5e3_Ca_0_1.mat');
lim3 = hstry.itr; pmem3 = hstry.pmem_track; pffil3 = hstry.pffil_track; pffir3 = hstry.pffir_track;
Nmem3 = hstry.Nmems; Nffil3 = hstry.Nffils; Nffir3 = hstry.Nffirs;
pmem3(:,:,2) = pmem3(:,:,2)-para.By; pffil3(:,:,2) = pffil3(:,:,2)-para.By; pffir3(:,:,2) = pffir3(:,:,2)-para.By;

load('wetting_20_10_gap_4_xl_2_gamma1_0_9_cb_5e3_Ca_0_1.mat');
lim4 = hstry.itr; pmem4 = hstry.pmem_track; pffil4 = hstry.pffil_track; pffir4 = hstry.pffir_track;
Nmem4 = hstry.Nmems; Nffil4 = hstry.Nffils; Nffir4 = hstry.Nffirs;
pmem4(:,:,2) = pmem4(:,:,2)-para.By; pffil4(:,:,2) = pffil4(:,:,2)-para.By; pffir4(:,:,2) = pffir4(:,:,2)-para.By;

end

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 0;

ind1 = floor(T/0.01)+1; %ind2 = floor(T/0.08)+1;

if ind1<lim1
    plot(pmem1(ind1,1:Nmem1(ind1),1),pmem1(ind1,1:Nmem1(ind1),2),'-.','color',c1,'LineWidth',3.0);
else
    plot(pmem1(lim1,1:Nmem1(lim1),1),pmem1(lim1,1:Nmem1(lim1),2),'-.','color',c1,'LineWidth',3.0);
end

hold on;
axis equal;
ylim([-0.7,0.7]);
xlim([0,para.Bx]);
set(gca,'YAxisLocation','right');

if ind1<lim2
    plot(pmem2(ind1,1:Nmem2(ind1),1),pmem2(ind1,1:Nmem2(ind1),2),'--','color',c2,'LineWidth',3.0);
else
    plot(pmem2(lim2,1:Nmem2(lim2),1),pmem2(lim2,1:Nmem2(lim2),2),'--','color',c2,'LineWidth',3.0);
end    
    
if ind1<lim3
    plot(pmem3(ind1,1:Nmem3(ind1),1),pmem3(ind1,1:Nmem3(ind1),2),':','color',c3,'LineWidth',3.0);
else
    plot(pmem3(lim3,1:Nmem3(lim3),1),pmem3(lim3,1:Nmem3(lim3),2),':','color',c3,'LineWidth',3.0);
end   

if ind1<lim4
    plot(pmem4(ind1,1:Nmem4(ind1),1),pmem4(ind1,1:Nmem4(ind1),2),'-','color',c4,'LineWidth',3.0);
else
    plot(pmem4(lim4,1:Nmem4(lim4),1),pmem4(lim4,1:Nmem4(lim4),2),'-','color',c4,'LineWidth',3.0);
end   
   
for plot_others = 1:1    
if ind1<lim1
    plot(pmem1(ind1,1:Nmem1(ind1),1),-pmem1(ind1,1:Nmem1(ind1),2),'-.','color',c1,'LineWidth',3.0);
    plot(pffil1(ind1,1:Nffil1(ind1),1),pffil1(ind1,1:Nffil1(ind1),2),'-.','color',c1,'LineWidth',3.0);
    plot(pffil1(ind1,1:Nffil1(ind1),1),-pffil1(ind1,1:Nffil1(ind1),2),'-.','color',c1,'LineWidth',3.0);
    plot(pffir1(ind1,1:Nffir1(ind1),1),pffir1(ind1,1:Nffir1(ind1),2),'-.','color',c1,'LineWidth',3.0);
    plot(pffir1(ind1,1:Nffir1(ind1),1),-pffir1(ind1,1:Nffir1(ind1),2),'-.','color',c1,'LineWidth',3.0);
else
    plot(pmem1(lim1,1:Nmem1(lim1),1),-pmem1(lim1,1:Nmem1(lim1),2),'-.','color',c1,'LineWidth',3.0);
    plot(pffil1(lim1,1:Nffil1(lim1),1),pffil1(lim1,1:Nffil1(lim1),2),'-.','color',c1,'LineWidth',3.0);
    plot(pffil1(lim1,1:Nffil1(lim1),1),-pffil1(lim1,1:Nffil1(lim1),2),'-.','color',c1,'LineWidth',3.0);
    plot(pffir1(lim1,1:Nffir1(lim1),1),pffir1(lim1,1:Nffir1(lim1),2),'-.','color',c1,'LineWidth',3.0); 
    plot(pffir1(lim1,1:Nffir1(lim1),1),-pffir1(lim1,1:Nffir1(lim1),2),'-.','color',c1,'LineWidth',3.0); 
end

if ind1<lim2
    plot(pmem2(ind1,1:Nmem2(ind1),1),-pmem2(ind1,1:Nmem2(ind1),2),'--','color',c2,'LineWidth',3.0);
    plot(pffil2(ind1,1:Nffil2(ind1),1),pffil2(ind1,1:Nffil2(ind1),2),'--','color',c2,'LineWidth',3.0);
    plot(pffil2(ind1,1:Nffil2(ind1),1),-pffil2(ind1,1:Nffil2(ind1),2),'--','color',c2,'LineWidth',3.0);
    plot(pffir2(ind1,1:Nffir2(ind1),1),pffir2(ind1,1:Nffir2(ind1),2),'--','color',c2,'LineWidth',3.0);
    plot(pffir2(ind1,1:Nffir2(ind1),1),-pffir2(ind1,1:Nffir2(ind1),2),'--','color',c2,'LineWidth',3.0);
else
    plot(pmem2(lim2,1:Nmem2(lim2),1),-pmem2(lim2,1:Nmem2(lim2),2),'--','color',c2,'LineWidth',3.0);
    plot(pffil2(lim2,1:Nffil2(lim2),1),pffil2(lim2,1:Nffil2(lim2),2),'--','color',c2,'LineWidth',3.0);
    plot(pffil2(lim2,1:Nffil2(lim2),1),-pffil2(lim2,1:Nffil2(lim2),2),'--','color',c2,'LineWidth',3.0);
    plot(pffir2(lim2,1:Nffir2(lim2),1),pffir2(lim2,1:Nffir2(lim2),2),'--','color',c2,'LineWidth',3.0); 
    plot(pffir2(lim2,1:Nffir2(lim2),1),-pffir2(lim2,1:Nffir2(lim2),2),'--','color',c2,'LineWidth',3.0); 
end

if ind1<lim3
    plot(pmem3(ind1,1:Nmem3(ind1),1),-pmem3(ind1,1:Nmem3(ind1),2),':','color',c3,'LineWidth',3.0);
    plot(pffil3(ind1,1:Nffil3(ind1),1),pffil3(ind1,1:Nffil3(ind1),2),':','color',c3,'LineWidth',3.0);
    plot(pffil3(ind1,1:Nffil3(ind1),1),-pffil3(ind1,1:Nffil3(ind1),2),':','color',c3,'LineWidth',3.0);
    plot(pffir3(ind1,1:Nffir3(ind1),1),pffir3(ind1,1:Nffir3(ind1),2),':','color',c3,'LineWidth',3.0);
    plot(pffir3(ind1,1:Nffir3(ind1),1),-pffir3(ind1,1:Nffir3(ind1),2),':','color',c3,'LineWidth',3.0);
else
    plot(pmem3(lim3,1:Nmem3(lim3),1),-pmem3(lim3,1:Nmem3(lim3),2),':','color',c3,'LineWidth',3.0);
    plot(pffil3(lim3,1:Nffil3(lim3),1),pffil3(lim3,1:Nffil3(lim3),2),':','color',c3,'LineWidth',3.0);
    plot(pffil3(lim3,1:Nffil3(lim3),1),-pffil3(lim3,1:Nffil3(lim3),2),':','color',c3,'LineWidth',3.0);
    plot(pffir3(lim3,1:Nffir3(lim3),1),pffir3(lim3,1:Nffir3(lim3),2),':','color',c3,'LineWidth',3.0); 
    plot(pffir3(lim3,1:Nffir3(lim3),1),-pffir3(lim3,1:Nffir3(lim3),2),':','color',c3,'LineWidth',3.0); 
end

if ind1<lim4
    plot(pmem4(ind1,1:Nmem4(ind1),1),-pmem4(ind1,1:Nmem4(ind1),2),'-','color',c4,'LineWidth',3.0);
    plot(pffil4(ind1,1:Nffil4(ind1),1),pffil4(ind1,1:Nffil4(ind1),2),'-','color',c4,'LineWidth',3.0);
    plot(pffil4(ind1,1:Nffil4(ind1),1),-pffil4(ind1,1:Nffil4(ind1),2),'-','color',c4,'LineWidth',3.0);
    plot(pffir4(ind1,1:Nffir4(ind1),1),pffir4(ind1,1:Nffir4(ind1),2),'-','color',c4,'LineWidth',3.0);
    plot(pffir4(ind1,1:Nffir4(ind1),1),-pffir4(ind1,1:Nffir4(ind1),2),'-','color',c4,'LineWidth',3.0);
else
    plot(pmem4(lim4,1:Nmem4(lim4),1),-pmem4(lim4,1:Nmem4(lim4),2),'-','color',c4,'LineWidth',3.0);
    plot(pffil4(lim4,1:Nffil4(lim4),1),pffil4(lim4,1:Nffil4(lim4),2),'-','color',c4,'LineWidth',3.0);
    plot(pffil4(lim4,1:Nffil4(lim4),1),-pffil4(lim4,1:Nffil4(lim4),2),'-','color',c4,'LineWidth',3.0);
    plot(pffir4(lim4,1:Nffir4(lim4),1),pffir4(lim4,1:Nffir4(lim4),2),'-','color',c4,'LineWidth',3.0); 
    plot(pffir4(lim4,1:Nffir4(lim4),1),-pffir4(lim4,1:Nffir4(lim4),2),'-','color',c4,'LineWidth',3.0); 
end 

end

plot([0 0],[-para.By, para.By],'-k','LineWidth',4.0); 


% hh = legend({'$\cos\theta_Y=0.7$','$\cos\theta_Y=0.5$',...
%                 '$\cos\theta_Y=0.3$','$\cos\theta_Y=0.1$'},'Interpreter','latex');
% set(hh,'FontName','Times New Roman','FontSize',32,'FontWeight','normal','Orientation','horizontal');

% title(['t = ',num2str((ind1-1)*0.04)]);
% xlabel('$x$','Interpreter','latex');
% ylabel('$y$','Interpreter','latex');
set(gca,'FontSize',32);

% plot([0 para.Bx], [-para.By -para.By],'-k','LineWidth',4.0);
% plot([0 para.Bx], [para.By para.By],'-k','LineWidth',4.0);


hold off;


