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
disp(['Time =',num2str(Time)]);

load('rectangle_dewetting_data32_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
plot(para.dt*(0:ind-1),-(para.areainit-hstry.droparea(1:ind))./para.areainit,...
    ':k','LineWidth',2.0);
hold on;

load('rectangle_dewetting_data64_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
plot(para.dt*(0:ind-1),-(para.areainit-hstry.droparea(1:ind))./para.areainit,...
    '--','LineWidth',2.0,'color',c2);

load('rectangle_dewetting_data128_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
plot(para.dt*(0:ind-1),-(para.areainit-hstry.droparea(1:ind))./para.areainit,...
    '-.','LineWidth',2.0,'color',c3);

load('rectangle_dewetting_data256_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
plot(para.dt*(0:ind-1),-(para.areainit-hstry.droparea(1:ind))./para.areainit,...
    '-','LineWidth',2.0,'color',c1);
disp(para.dt);

leg = legend('$\mathcal{Q}=32$','$\mathcal{Q}=64$', '$\mathcal{Q}=128$', '$\mathcal{Q}=256$');
set(leg,'Interpreter','latex'); 
xlabel('$t$','Interpreter','latex');
ylabel('$\Delta A(t)$','Interpreter','latex');
set(gca,'FontSize',32);
hold off;

