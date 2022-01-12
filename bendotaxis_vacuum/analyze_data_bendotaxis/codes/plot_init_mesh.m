close all;
clear;
clc;
addpath(genpath('.'));

load('init_vacuum_testconv_prewetting_20_3_0_5_gap_1_xl_1_gamma1_0_5.mat');

elem = geom.elem; node = geom.node;
pmem = sln.pmem; pffil = sln.pffil; pffir = sln.pffir;
Nmem = para.Nmem; Nffil = para.Nffil; Nffir = para.Nffir;

trimesh(elem,node(:,1),node(:,2),zeros(size(node,1),1),'facecolor','w','edgecolor','k','LineWidth',1.0);
hold on;
plot([0, para.Bx],[para.By para.By],'-k','LineWidth',1.5);
plot(pmem(1:Nmem,1),pmem(1:Nmem,2),'-bs','LineWidth',1.5,'MarkerSize',3,'MarkerFaceColor','b');
plot(pffil(1:Nffil,1),pffil(1:Nffil,2),'-ro','LineWidth',1.5,'MarkerSize',2,'MarkerFaceColor','r');
plot(pffir(1:Nffir,1),pffir(1:Nffir,2),'-ro','LineWidth',1.5,'MarkerSize',2,'MarkerFaceColor','r');
view(2);
axis equal;
% axis off;
hold off;