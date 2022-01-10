close all;
clear;
clc;
addpath(genpath('.'));
% load('initdata_adapt_256_cl.mat');
load('initdata_adapt_64_cl_rectangle.mat');

% elem = geom.elem; node = geom.node;
pmem = sln.pmem; pffi = sln.pffi;
Nmem = para.Nmem; Nffi = para.Nffi;
Bx = para.Bx; By = para.By;

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end   

plot(pmem(1:Nmem,1),pmem(1:Nmem,2),'-','LineWidth',3.0,'color',c3);
hold on;
plot(pffi(1:Nffi,1),pffi(1:Nffi,2),'-r','LineWidth',3.0);
plot([-Bx -Bx],[pmem(1,2) By],'-k','LineWidth',1.5);
plot([Bx Bx],[pmem(Nmem,2) By],'-k','LineWidth',1.5);
plot([-Bx Bx],[By By],'-k','LineWidth',1.5);
axis equal;
% xlim([-1.03,1.03]);
% ylim([-0.1,1]);
xlim([-1,1]);
ylim([-0.15,1]);
set(gca,'FontSize',32);
hold off;