clear;
clc;
addpath(genpath('.'));
dbstop if error


load('dewetting_20_10_gap_4_xl_2_gamma2_0_3_cb_5e3_Ca_0_1.mat');

dt = para.dt; cl2mems = hstry.cl2mem;
pmems = hstry.pmem_track; nus = hstry.nu_track;

for i = 2:100:hstry.itr
    Nmem = hstry.Nmems(i); cl2mem1 = cl2mems(i,1); cl2mem2 = cl2mems(i,2); Nm = hstry.Nms(i);
    plot(pmems(i-1,cl2mem1:cl2mem2,1),nus(i,1:Nm),'-ob');
    hold on;
    plot([pmems(i-1,1,1) pmems(i-1,cl2mem1,1)],[nus(i,Nm+1) nus(i,Nm+1)],'-r');
    plot([pmems(i-1,cl2mem2,1) pmems(i-1,Nmem,1)], [0 0], '-r');
    hh = legend(['t = ',num2str((i-1)*dt)],...
           ['left nu jump: ', num2str(nus(i,1)-nus(i,Nm+1))],...
           ['right nu jump: ', num2str(nus(i,Nm))]);
    set(hh,'FontName','Times New Roman','FontSize',18,'FontWeight','normal');   
%     axis equal;
%     axis tight;
    hold off;
    pause(0.001);
end
