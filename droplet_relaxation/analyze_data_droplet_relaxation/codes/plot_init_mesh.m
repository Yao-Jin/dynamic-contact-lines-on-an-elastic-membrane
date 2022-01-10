close all;
clear;
clc;
addpath(genpath('.'));
load('initdata_adapt_64_cl_rectangle.mat');

elem = geom.elem; node = geom.node;
pmem = sln.pmem; pffi = sln.pffi;
Nmem = para.Nmem; Nffi = para.Nffi;

trimesh(elem,node(:,1),node(:,2),zeros(size(node,1),1),'facecolor','w','edgecolor','k','LineWidth',1.0);
hold on;
plot(pmem(1:Nmem,1),pmem(1:Nmem,2),'-bs','LineWidth',1.5,'MarkerSize',5,'MarkerFaceColor','b');
plot(pffi(1:Nffi,1),pffi(1:Nffi,2),'-ro','LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor','r');
view(2);
axis equal;
axis off;
hold off;