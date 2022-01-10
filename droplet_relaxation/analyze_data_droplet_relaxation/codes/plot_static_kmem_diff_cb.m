% close all;
clear;
% clc;
addpath(genpath('.'));
dbstop if error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time = 2; lw = 2;

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end 

load('rectangle_wetting_data256_cb1e_3_gamma1_0_5_reprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
p1 = plot(hstry.pmem_track(ind-1,1:hstry.Nmems(ind-1),1),...
     hstry.kmem_track(ind,1:hstry.Nmems(ind)),'-.o','LineWidth',lw,'color',c1);
hold on;
 
load('rectangle_dewetting_data256_cb1e_3_gamma2_0_5_reprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
p2 = plot(hstry.pmem_track(ind-1,1:hstry.Nmems(ind-1),1),...
     hstry.kmem_track(ind,1:hstry.Nmems(ind)),'-o','LineWidth',lw,'color',c2);

leg = legend([p1,p2],{'$\theta_Y=\pi/3,\ c_b=0.001$','$\theta_Y=2\pi/3,\ c_b=0.001$'});
set(leg,'Interpreter','latex'); 
xlabel('$x$','Interpreter','latex');
ylabel('$\kappa$','Interpreter','latex');
set(gca,'FontSize',32);
% axis equal;
ylim([-17,2.5]); 
xlim([-1,1]);
hold off;
