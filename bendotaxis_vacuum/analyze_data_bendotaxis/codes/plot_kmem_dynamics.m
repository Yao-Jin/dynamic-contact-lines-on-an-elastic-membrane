clear;
clc;
addpath(genpath('.'));
dbstop if error

load('wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1.mat');


dt = para.dt; cl2mems = hstry.cl2mem;
pmems = hstry.pmem_track; kmems = hstry.kmem_track;

for i = 2:16:hstry.itr
    Nmem = hstry.Nmems(i); cl2mem1 = cl2mems(i,1); cl2mem2 = cl2mems(i,2); Nm = hstry.Nms(i);
    plot(pmems(i-1,1:Nmem,1),kmems(i,1:Nmem),'-ob');
    hold on;
    hh = legend(['t = ',num2str((i-1)*dt)]);
    set(hh,'FontName','Times New Roman','FontSize',18,'FontWeight','normal');   
%     axis equal;
%     axis tight;
    hold off;
    pause(0.001);
end
