clear;
clc;
addpath(genpath('.'));
dbstop if error

% load('rectangle_dewetting_data128_cb1e_4_gamma2_0_5_reprod_Ca_0_2.mat');
load('rectangle_dewetting_data256_cb1e_3_gamma2_0_5_reprod_Ca_0_2.mat');

dt = para.dt;
pmems = hstry.pmem_track; pffis = hstry.pffi_track;

for i = 1:100:hstry.itr
    Nmem = hstry.Nmems(i); Nffi = hstry.Nffis(i);
    plot(pmems(i,1:Nmem,1),pmems(i,1:Nmem,2),'-ok');
    hold on;
    plot(pffis(i,1:Nffi,1),pffis(i,1:Nffi,2),'-sr');
    hh = legend(['t = ',num2str((i-1)*dt)]);
    set(hh,'FontName','Times New Roman','FontSize',18,'FontWeight','normal');
    axis equal;
    hold off;
    pause(0.01);
end
