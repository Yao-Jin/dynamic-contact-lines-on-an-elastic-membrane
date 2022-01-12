clear;
clc;
addpath(genpath('.'));
dbstop if error

load('wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1.mat');
% load('dewetting_20_10_gap_4_xl_2_gamma2_0_3_cb_5e3_Ca_0_1.mat');

dt = para.dt; Nffil = para.Nffil; Nffir = para.Nffir;
pmems = hstry.pmem_track; pffils = hstry.pffil_track; pffirs = hstry.pffir_track;

for i = 1:100:hstry.itr
    Nmem = hstry.Nmems(i); Nffil = hstry.Nffils(i); Nffir = hstry.Nffirs(i); 
    plot(pmems(i,1:Nmem,1),pmems(i,1:Nmem,2),'-ob');
    hold on;
    plot([0 para.Bx], [para.By para.By],'-k');
    plot(pffils(i,1:Nffil,1),pffils(i,1:Nffil,2),'-sr');
    plot(pffirs(i,1:Nffir,1),pffirs(i,1:Nffir,2),'-sr');
    hh = legend(['t = ',num2str((i-1)*dt)]);
    set(hh,'FontName','Times New Roman','FontSize',18,'FontWeight','normal');
    axis equal;
    hold off;
    pause(0.001);
end
