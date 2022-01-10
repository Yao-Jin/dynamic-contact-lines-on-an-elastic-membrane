% close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end 

% load('record_curvature_rect_wetting_data128.mat');
load('record_curvature_rect_dewetting_data128.mat');

ind = floor(0.05/para.dt)+1; Nmem = hstry.Nmems(ind);
cl2mem1 = hstry.cl2mem(ind,1); cl2mem2 = hstry.cl2mem(ind,2);
pmem = reshape(hstry.pmem_track(ind-1,1:2*Nmem-1,1:2),2*Nmem-1,2);
nu = reshape(hstry.nu_track(ind,1:2*Nmem+1),2*Nmem+1,1);
    pmem2 = zeros(2*Nmem-1,2);
    pmem2(1:2:2*Nmem-1,1:2) = pmem(1:Nmem,1:2);
    pmem2(2:2:2*Nmem-2,1:2) = pmem(Nmem+1:2*Nmem-1,1:2); 
    nu2 = zeros(2*Nmem+1,1);
    nu2(1:2:2*cl2mem1-1) = nu(1:cl2mem1);
    nu2(2:2:2*cl2mem1-2) = nu(Nmem+3:Nmem+2+cl2mem1-1);
    nu2(2*cl2mem1-1:2:2*cl2mem2-1) = nu(cl2mem1:cl2mem2);
    nu2(2*cl2mem1:2:2*cl2mem2-2)= nu(Nmem+2+cl2mem1:Nmem+2+cl2mem2-1);
    nu2(2*cl2mem2-1:2:2*Nmem-1) = nu(cl2mem2:Nmem);
    nu2(2*cl2mem2:2:2*Nmem-2) = nu(Nmem+2+cl2mem2:Nmem+2+Nmem-1);
    
p1 = plot(pmem2(1:2*cl2mem1-1,1),nu2(1:2*cl2mem1-1),'-','LineWidth',2.0,'color',c1);
hold on;
plot(pmem2(2*cl2mem1-1:2*cl2mem2-1,1),[nu(Nmem+1); nu2(2*cl2mem1:2*cl2mem2-2); nu(Nmem+2)],...
     '-','LineWidth',2.0,'color',c1);
plot(pmem2(2*cl2mem2-1:2*Nmem-2,1),nu2(2*cl2mem2-1:2*Nmem-2),'-','LineWidth',2.0,'color',c1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind = floor(0.2/para.dt)+1;
cl2mem1 = hstry.cl2mem(ind,1); cl2mem2 = hstry.cl2mem(ind,2);
pmem = reshape(hstry.pmem_track(ind-1,1:2*Nmem-1,1:2),2*Nmem-1,2);
nu = reshape(hstry.nu_track(ind,1:2*Nmem+1),2*Nmem+1,1);
    pmem2 = zeros(2*Nmem-1,2);
    pmem2(1:2:2*Nmem-1,1:2) = pmem(1:Nmem,1:2);
    pmem2(2:2:2*Nmem-2,1:2) = pmem(Nmem+1:2*Nmem-1,1:2); 
    nu2 = zeros(2*Nmem+1,1);
    nu2(1:2:2*cl2mem1-1) = nu(1:cl2mem1);
    nu2(2:2:2*cl2mem1-2) = nu(Nmem+3:Nmem+2+cl2mem1-1);
    nu2(2*cl2mem1-1:2:2*cl2mem2-1) = nu(cl2mem1:cl2mem2);
    nu2(2*cl2mem1:2:2*cl2mem2-2)= nu(Nmem+2+cl2mem1:Nmem+2+cl2mem2-1);
    nu2(2*cl2mem2-1:2:2*Nmem-1) = nu(cl2mem2:Nmem);
    nu2(2*cl2mem2:2:2*Nmem-2) = nu(Nmem+2+cl2mem2:Nmem+2+Nmem-1);
    
p2 = plot(pmem2(1:2*cl2mem1-1,1),nu2(1:2*cl2mem1-1),'--','LineWidth',2.0,'color',c2);
plot(pmem2(2*cl2mem1-1:2*cl2mem2-1,1),[nu(Nmem+1); nu2(2*cl2mem1:2*cl2mem2-2); nu(Nmem+2)],...
     '--','LineWidth',2.0,'color',c2);
plot(pmem2(2*cl2mem2-1:2*Nmem-2,1),nu2(2*cl2mem2-1:2*Nmem-2),'--','LineWidth',2.0,'color',c2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind = floor(0.5/para.dt)+1;
cl2mem1 = hstry.cl2mem(ind,1); cl2mem2 = hstry.cl2mem(ind,2);
pmem = reshape(hstry.pmem_track(ind-1,1:2*Nmem-1,1:2),2*Nmem-1,2);
nu = reshape(hstry.nu_track(ind,1:2*Nmem+1),2*Nmem+1,1);
    pmem2 = zeros(2*Nmem-1,2);
    pmem2(1:2:2*Nmem-1,1:2) = pmem(1:Nmem,1:2);
    pmem2(2:2:2*Nmem-2,1:2) = pmem(Nmem+1:2*Nmem-1,1:2); 
    nu2 = zeros(2*Nmem+1,1);
    nu2(1:2:2*cl2mem1-1) = nu(1:cl2mem1);
    nu2(2:2:2*cl2mem1-2) = nu(Nmem+3:Nmem+2+cl2mem1-1);
    nu2(2*cl2mem1-1:2:2*cl2mem2-1) = nu(cl2mem1:cl2mem2);
    nu2(2*cl2mem1:2:2*cl2mem2-2)= nu(Nmem+2+cl2mem1:Nmem+2+cl2mem2-1);
    nu2(2*cl2mem2-1:2:2*Nmem-1) = nu(cl2mem2:Nmem);
    nu2(2*cl2mem2:2:2*Nmem-2) = nu(Nmem+2+cl2mem2:Nmem+2+Nmem-1);
    
p3 = plot(pmem2(1:2*cl2mem1-1,1),nu2(1:2*cl2mem1-1),'-.','LineWidth',2.0,'color',c3);
plot(pmem2(2*cl2mem1-1:2*cl2mem2-1,1),[nu(Nmem+1); nu2(2*cl2mem1:2*cl2mem2-2); nu(Nmem+2)],...
     '-.','LineWidth',2.0,'color',c3);
plot(pmem2(2*cl2mem2-1:2*Nmem-2,1),nu2(2*cl2mem2-1:2*Nmem-2),'-.','LineWidth',2.0,'color',c3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind = floor(1/para.dt)+1;
cl2mem1 = hstry.cl2mem(ind,1); cl2mem2 = hstry.cl2mem(ind,2);
pmem = reshape(hstry.pmem_track(ind-1,1:2*Nmem-1,1:2),2*Nmem-1,2);
nu = reshape(hstry.nu_track(ind,1:2*Nmem+1),2*Nmem+1,1);
    pmem2 = zeros(2*Nmem-1,2);
    pmem2(1:2:2*Nmem-1,1:2) = pmem(1:Nmem,1:2);
    pmem2(2:2:2*Nmem-2,1:2) = pmem(Nmem+1:2*Nmem-1,1:2); 
    nu2 = zeros(2*Nmem+1,1);
    nu2(1:2:2*cl2mem1-1) = nu(1:cl2mem1);
    nu2(2:2:2*cl2mem1-2) = nu(Nmem+3:Nmem+2+cl2mem1-1);
    nu2(2*cl2mem1-1:2:2*cl2mem2-1) = nu(cl2mem1:cl2mem2);
    nu2(2*cl2mem1:2:2*cl2mem2-2)= nu(Nmem+2+cl2mem1:Nmem+2+cl2mem2-1);
    nu2(2*cl2mem2-1:2:2*Nmem-1) = nu(cl2mem2:Nmem);
    nu2(2*cl2mem2:2:2*Nmem-2) = nu(Nmem+2+cl2mem2:Nmem+2+Nmem-1);
    
p4 = plot(pmem2(1:2*cl2mem1-1,1),nu2(1:2*cl2mem1-1),'-','LineWidth',2.0,'color',c4);
plot(pmem2(2*cl2mem1-1:2*cl2mem2-1,1),[nu(Nmem+1); nu2(2*cl2mem1:2*cl2mem2-2); nu(Nmem+2)],...
     '-','LineWidth',2.0,'color',c4);
plot(pmem2(2*cl2mem2-1:2*Nmem-2,1),nu2(2*cl2mem2-1:2*Nmem-2),'-','LineWidth',2.0,'color',c4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind = floor(2/para.dt)+1;
cl2mem1 = hstry.cl2mem(ind,1); cl2mem2 = hstry.cl2mem(ind,2);
pmem = reshape(hstry.pmem_track(ind-1,1:2*Nmem-1,1:2),2*Nmem-1,2);
nu = reshape(hstry.nu_track(ind,1:2*Nmem+1),2*Nmem+1,1);
    pmem2 = zeros(2*Nmem-1,2);
    pmem2(1:2:2*Nmem-1,1:2) = pmem(1:Nmem,1:2);
    pmem2(2:2:2*Nmem-2,1:2) = pmem(Nmem+1:2*Nmem-1,1:2); 
    nu2 = zeros(2*Nmem+1,1);
    nu2(1:2:2*cl2mem1-1) = nu(1:cl2mem1);
    nu2(2:2:2*cl2mem1-2) = nu(Nmem+3:Nmem+2+cl2mem1-1);
    nu2(2*cl2mem1-1:2:2*cl2mem2-1) = nu(cl2mem1:cl2mem2);
    nu2(2*cl2mem1:2:2*cl2mem2-2)= nu(Nmem+2+cl2mem1:Nmem+2+cl2mem2-1);
    nu2(2*cl2mem2-1:2:2*Nmem-1) = nu(cl2mem2:Nmem);
    nu2(2*cl2mem2:2:2*Nmem-2) = nu(Nmem+2+cl2mem2:Nmem+2+Nmem-1);
    
p5 = plot(pmem2(1:2*cl2mem1-1,1),nu2(1:2*cl2mem1-1),':','LineWidth',2.0,'color',c5);
plot(pmem2(2*cl2mem1-1:2*cl2mem2-1,1),[nu(Nmem+1); nu2(2*cl2mem1:2*cl2mem2-2); nu(Nmem+2)],...
     ':','LineWidth',2.0,'color',c5);
plot(pmem2(2*cl2mem2-1:2*Nmem-2,1),nu2(2*cl2mem2-1:2*Nmem-2),':','LineWidth',2.0,'color',c5);

leg = legend([p1,p2,p3,p4,p5],{'$t=0.05$','$t=0.2$','$t=0.5$','$t=1$','$t=2$'});
set(leg,'Interpreter','latex','Location','southeast'); 
xlabel('$x$','Interpreter','latex');
ylabel('$\nu$','Interpreter','latex');
set(gca,'FontSize',32);
% axis equal;
% ylim([-2.5,0]);
% xlim([-1,1]);
hold off;
