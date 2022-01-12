% close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error


% cos(thetaY) = (gamma2-gamma1)/gamma3; gamma3 = 1;
% wetting: gamma2 = 1;
% dewetting: gamma1 = 1;

load('wetting_20_10_gap_4_xl_2_gamma1_0_5_cb_1e3_Ca_0_1.mat');
lim = hstry.itr; pmem = hstry.pmem_track; pffil = hstry.pffil_track; pffir = hstry.pffir_track;
Nmem = hstry.Nmems; Nffil = hstry.Nffils; Nffir = hstry.Nffirs;
pmem(:,:,2) = pmem(:,:,2)-para.By; pffil(:,:,2) = pffil(:,:,2)-para.By; pffir(:,:,2) = pffir(:,:,2)-para.By;

c = zeros(8,3);
for definecolor = 1:1
    c(1,1:3) = [228 26 28]/256; c(2,1:3) = [55 126 184]/256; c(3,1:3) = [77 175 74]/256; c(4,1:3) = [152 78 163]/256;
    c(5,1:3) = [255 127 0]/256; c(6,1:3) = [255 217 47]/256; c(7,1:3) = [166 86 40]/256; c(8,1:3) = [246 112 136]/256;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot([0 0],[-para.By, para.By],'-k','LineWidth',4.0); 
hold on;
axis equal;
ylim([-0.7,0.7]);
xlim([0,para.Bx]);
set(gca,'FontSize',32);
set(gca,'YAxisLocation','right');

Times = [0 4 8 11]/10; linetypes = [":"; "-."; "--"; "-"];
for i=1:4
T = Times(i); color_type = c(i,1:3); linetype = linetypes(i);  
ind1 = floor(T/para.dt)+1;

if ind1<lim
    plot(pmem(ind1,1:Nmem(ind1),1),pmem(ind1,1:Nmem(ind1),2),linetype,'color',color_type,'LineWidth',3.0);
else
    plot(pmem(lim,1:Nmem(lim),1),pmem(lim,1:Nmem(lim),2),linetype,'color',color_type,'LineWidth',3.0);
end    
   
if ind1<lim
    plot(pmem(ind1,1:Nmem(ind1),1),-pmem(ind1,1:Nmem(ind1),2),linetype,'color',color_type,'LineWidth',3.0);
    plot(pffil(ind1,1:Nffil(ind1),1),pffil(ind1,1:Nffil(ind1),2),linetype,'color',color_type,'LineWidth',3.0);
    plot(pffil(ind1,1:Nffil(ind1),1),-pffil(ind1,1:Nffil(ind1),2),linetype,'color',color_type,'LineWidth',3.0);
    plot(pffir(ind1,1:Nffir(ind1),1),pffir(ind1,1:Nffir(ind1),2),linetype,'color',color_type,'LineWidth',3.0);
    plot(pffir(ind1,1:Nffir(ind1),1),-pffir(ind1,1:Nffir(ind1),2),linetype,'color',color_type,'LineWidth',3.0);
else
    plot(pmem(lim,1:Nmem(lim),1),-pmem(lim,1:Nmem(lim),2),linetype,'color',color_type,'LineWidth',3.0);
    plot(pffil(lim,1:Nffil(lim),1),pffil(lim,1:Nffil(lim),2),linetype,'color',color_type,'LineWidth',3.0);
    plot(pffil(lim,1:Nffil(lim),1),-pffil(lim,1:Nffil(lim),2),linetype,'color',color_type,'LineWidth',3.0);
    plot(pffir(lim,1:Nffir(lim),1),pffir(lim,1:Nffir(lim),2),linetype,'color',color_type,'LineWidth',3.0); 
    plot(pffir(lim,1:Nffir(lim),1),-pffir(lim,1:Nffir(lim),2),linetype,'color',color_type,'LineWidth',3.0); 
end

end

% hh = legend({'$t=0$','$t=4$','$t=8$','$t=11$'},'Interpreter','latex');
% set(hh,'FontSize',32,'Orientation','horizontal');

hold off;


