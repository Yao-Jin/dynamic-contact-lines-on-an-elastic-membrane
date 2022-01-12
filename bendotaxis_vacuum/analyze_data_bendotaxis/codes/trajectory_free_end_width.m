clear;
clc;
addpath(genpath('.'));
dbstop if error

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end

load('wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1.mat');
dt = para.dt; %By = para.By;
pmems = hstry.pmem_track;
free_end1 = zeros(hstry.itr,1);
for i = 1:1:hstry.itr
    Nmem = hstry.Nmems(i);
    free_end1(i) = pmems(i,Nmem,2);
end

load('dewetting_20_10_gap_4_xl_2_gamma2_0_3_cb_5e3_Ca_0_1.mat');
T =50;
ind = round(T/para.dt)+1;
pmems = hstry.pmem_track;
free_end2 = zeros(ind,1);
for i = 1:1:ind
    Nmem = hstry.Nmems(i);
    free_end2(i) = pmems(i,Nmem,2);
end

T =16;
ind1 = round(T/para.dt)+1;
T =16;
ind2 = round(T/para.dt)+1;

p1 = plot((0:ind1-1)*dt,1-2*free_end1(1:ind1),'-r','LineWidth',2.0);
hold on;
p2 = plot((0:ind2-1)*dt,1-2*free_end2(1:ind2),'-b','LineWidth',2.0);
% plot([0.4 0.4],[0.25 -0.1],'-k','LineWidth',1.0);
plot([0.3 0.3],[1.2 0.4],'-k','LineWidth',1.0);
hh = legend([p1,p2],{'$\cos\theta_Y=0.7$','$\cos\theta_Y=-0.7$'});
set(hh,'Interpreter','latex','FontName','Times New Roman','FontSize',32,'FontWeight','normal','Location','northeast');
xlim([0,16]);
ylim([0.4,1.5]);
xlabel('$t$','Interpreter','latex');
ylabel('$w(t)$','Interpreter','latex');
set(gca,'FontSize',32);
% set(gca,'XAxisLocation','origin');
hold off;

axes('position',[0.18,0.15,0.25,0.25]);
load('record_wetting_gap_1_xl_1_gamma1_0_5_cb_5e1_Ca_1_t_3.mat');
start = 1; Nmem = para.Nmem; Nffil = para.Nffil; Nffir = para.Nffir;
pmem = sln.pmemold; pffil = sln.pffilold; pffir = sln.pffirold; By = para.By;
lw = 1.5;
plot(pmem(start:Nmem,1),pmem(start:Nmem,2)-By,'-','LineWidth',lw,'color',c2);
hold on;
plot(pffil(1:Nffil,1),pffil(1:Nffil,2)-By,'-','LineWidth',lw,'color',c1);
plot(pffir(1:Nffir,1),pffir(1:Nffir,2)-By,'-','LineWidth',lw,'color',c1);
plot(pmem(start:Nmem,1),By-pmem(start:Nmem,2),'-','LineWidth',lw,'color',c2);
plot(pffil(1:Nffil,1),By-pffil(1:Nffil,2),'-','LineWidth',lw,'color',c1);
plot(pffir(1:Nffir,1),By-pffir(1:Nffir,2),'-','LineWidth',lw,'color',c1);
plot([pmem(start,1) pmem(start,1)],[pmem(start,2)-By 0],'-k','LineWidth',3.0);
plot([pmem(start,1) pmem(start,1)],[pmem(start,2)+By 0],'-k','LineWidth',3.0);
axis equal;
axis off;
hold off;

axes('position',[0.18,0.65,0.25,0.25]);
load('record_dewetting_gap_2_xl_0_5_gamma2_0_1_cb_1e1_Ca_1_t_14.mat');
start = 1; Nmem = para.Nmem; Nffil = para.Nffil; Nffir = para.Nffir;
pmem = sln.pmemold; pffil = sln.pffilold; pffir = sln.pffirold; By = para.By;
lw = 1.5;
plot(pmem(start:Nmem,1),pmem(start:Nmem,2)-By,'-','LineWidth',lw,'color',c2);
hold on;
plot(pffil(1:Nffil,1),pffil(1:Nffil,2)-By,'-','LineWidth',lw,'color',c1);
plot(pffir(1:Nffir,1),pffir(1:Nffir,2)-By,'-','LineWidth',lw,'color',c1);
plot(pmem(start:Nmem,1),By-pmem(start:Nmem,2),'-','LineWidth',lw,'color',c2);
plot(pffil(1:Nffil,1),By-pffil(1:Nffil,2),'-','LineWidth',lw,'color',c1);
plot(pffir(1:Nffir,1),By-pffir(1:Nffir,2),'-','LineWidth',lw,'color',c1);
plot([pmem(start,1) pmem(start,1)],[pmem(start,2)-By 0],'-k','LineWidth',3.0);
plot([pmem(start,1) pmem(start,1)],[pmem(start,2)+By 0],'-k','LineWidth',3.0);
axis equal;
axis off;
hold off;
