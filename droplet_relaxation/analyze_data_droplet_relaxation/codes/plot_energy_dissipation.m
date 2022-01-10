close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error
 
load('rectangle_wetting_data128_cb1e_1_gamma1_0_5_reprod_Ca_0_2.mat');

energy_wetting = zeros(hstry.itr-1,1); willmore_wetting = zeros(hstry.itr-1,1);
for ind = 1:hstry.itr
    Nmem = hstry.Nmems(ind); Nffi = para.Nffi;
    pmem = reshape(hstry.pmem_track(ind,1:Nmem,1:2),Nmem,2); 
    pffi = reshape(hstry.pffi_track(ind,1:Nffi,1:2),Nffi,2);
    kmem_nodal = getcurvature(pmem); kmem_nodal(1) = 0; kmem_nodal(Nmem) = 0;
    % trapezoidal rule
    willmore = 0;
    x = hstry.pmem_track(ind,1:Nmem,1); y = hstry.pmem_track(ind,1:Nmem,2);
    for i = 1:Nmem-1
        loclen = sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
        willmore = willmore + (kmem_nodal(i)^2+kmem_nodal(i+1)^2)/2*loclen; 
    end
%     cl2mem1 = hstry.cl2mem(ind,1); cl2mem2 = hstry.cl2mem(ind,2);
    for i = 1:Nmem
        if abs(pmem(i,1)-pffi(1,1))<1e-8
            cl2mem1 = i;
        end
        if abs(pmem(i,1)-pffi(Nffi,1))<1e-8
            cl2mem2 = i;
        end
    end

    energy_wetting(ind) = para.cb/2*willmore...
                  + para.gamma1*getlength(pmem(cl2mem1:cl2mem2,1:2))...
                  + para.gamma2*getlength(pmem(1:cl2mem1,1:2))...
                  + para.gamma2*getlength(pmem(cl2mem2:Nmem,1:2))...
                  + para.gamma3*getlength(pffi(1:Nffi,1:2));
    willmore_wetting(ind) = para.cb/2*willmore;
end

load('rectangle_dewetting_data128_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');

energy_dewetting = zeros(hstry.itr-1,1);  willmore_dewetting = zeros(hstry.itr-1,1);
for ind = 1:hstry.itr
    Nmem = hstry.Nmems(ind); Nffi = para.Nffi;
    pmem = reshape(hstry.pmem_track(ind,1:Nmem,1:2),Nmem,2); 
    pffi = reshape(hstry.pffi_track(ind,1:Nffi,1:2),Nffi,2);
    kmem_nodal = getcurvature(pmem); kmem_nodal(1) = 0; kmem_nodal(Nmem) = 0;
    % trapezoidal rule
    willmore = 0;
    x = hstry.pmem_track(ind,1:Nmem,1); y = hstry.pmem_track(ind,1:Nmem,2);
    for i = 1:Nmem-1
        loclen = sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
        willmore = willmore + (kmem_nodal(i)^2+kmem_nodal(i+1)^2)/2*loclen; 
    end
%     cl2mem1 = hstry.cl2mem(ind,1); cl2mem2 = hstry.cl2mem(ind,2);
    for i = 1:Nmem
        if abs(pmem(i,1)-pffi(1,1))<1e-8
            cl2mem1 = i;
        end
        if abs(pmem(i,1)-pffi(Nffi,1))<1e-8
            cl2mem2 = i;
        end
    end

    energy_dewetting(ind) = para.cb/2*willmore...
                  + para.gamma1*getlength(pmem(cl2mem1:cl2mem2,1:2))...
                  + para.gamma2*getlength(pmem(1:cl2mem1,1:2))...
                  + para.gamma2*getlength(pmem(cl2mem2:Nmem,1:2))...
                  + para.gamma3*getlength(pffi(1:Nffi,1:2));
    willmore_dewetting(ind) = para.cb/2*willmore;           
end

plot(para.dt*(0:hstry.itr-1),energy_wetting./energy_wetting(1),'-r','LineWidth',2.0);
hold on;
plot(para.dt*(0:hstry.itr-1),energy_dewetting./energy_dewetting(1),'--b','LineWidth',2.0);
legend('\theta_Y=\pi/3', '\theta_Y=2\pi/3');
xlim([0,2]);
xlabel('$t$','Interpreter','latex');
ylabel('$\mathcal{E}(t)/\mathcal{E}(0)$','Interpreter','latex');
set(gca,'FontSize',32);
hold off;

plot(para.dt*(0:hstry.itr-1),willmore_wetting,'-r','LineWidth',2.0);
hold on;
plot(para.dt*(0:hstry.itr-1),willmore_dewetting,'--b','LineWidth',2.0);
legend('\theta_Y=\pi/3', '\theta_Y=2\pi/3');
xlim([0,2]);
xlabel('$t$','Interpreter','latex');
ylabel('$\mathcal{E}_w(t)$','Interpreter','latex');
set(gca,'FontSize',32);
hold off;