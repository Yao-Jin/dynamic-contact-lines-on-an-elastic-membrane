close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end 

Time = 2;

load('rectangle_dewetting_data32_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
plot(para.dt*(0:ind-1),hstry.dynamic_angle_l(1:ind),':k','LineWidth',2.0);
hstry.dynamic_angle_l(ind)-para.thetaY
hold on;

load('rectangle_dewetting_data64_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
plot(para.dt*(0:ind-1),hstry.dynamic_angle_l(1:ind),'--','LineWidth',2.0,'color',c2);
hstry.dynamic_angle_l(ind)-para.thetaY

load('rectangle_dewetting_data128_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
plot(para.dt*(0:ind-1),hstry.dynamic_angle_l(1:ind),'-.','LineWidth',2.0,'color',c3);
hstry.dynamic_angle_l(ind)-para.thetaY

load('rectangle_dewetting_data256_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
plot(para.dt*(0:ind-1),hstry.dynamic_angle_l(1:ind),'-','LineWidth',2.0,'color',c1);
hstry.dynamic_angle_l(ind)-para.thetaY

plot([0 Time],[para.thetaY para.thetaY],'-','LineWidth',2.0,'color',c4);
leg = legend('$\mathcal{Q}=32$','$\mathcal{Q}=64$', '$\mathcal{Q}=128$', '$\mathcal{Q}=256$', '$\theta_Y$');
set(leg,'Interpreter','latex'); 
xlabel('$t$','Interpreter','latex');
ylabel('$\theta_d$','Interpreter','latex');
set(gca,'FontSize',32);
ylim([1.55,2.1]);
% axis tight;
% axis equal;
hold off;
