close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end 

Time = 3;
disp(['Time =',num2str(Time)]);

% load('test.mat');
% load('testconv4_10_cb_1e2.mat');
load('testconv_10_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');
ind = floor(Time/para.dt)+1;
plot(para.dt*(0:ind-1),-(para.areainit-hstry.droparea(1:ind))./para.areainit,...
    ':k','LineWidth',2.0);
% disp(para.dt);
hold on;

% load('testconv4_20_cb_1e2.mat');
load('testconv_20_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');
ind = floor(Time/para.dt)+1;
plot(para.dt*(0:ind-1),-(para.areainit-hstry.droparea(1:ind))./para.areainit,...
    '--','LineWidth',2.0,'color',c2);
% disp(para.dt);

% load('testconv4_40_cb_1e2.mat');
load('testconv_40_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');
ind = floor(Time/para.dt)+1;
plot(para.dt*(0:ind-1),-(para.areainit-hstry.droparea(1:ind))./para.areainit,...
    '-.','LineWidth',2.0,'color',c3);
% disp(para.dt);

% load('testconv4_80_cb_1e2.mat');
load('testconv_80_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');
ind = floor(Time/para.dt)+1;
plot(para.dt*(0:ind-1),-(para.areainit-hstry.droparea(1:ind))./para.areainit,...
    '-','LineWidth',2.0,'color',c1);
disp(para.dt);

leg = legend('$\mathcal{R}=10$','$\mathcal{R}=20$', '$\mathcal{R}=40$', '$\mathcal{R}=80$');
set(leg,'Interpreter','latex'); 
xlabel('$t$','Interpreter','latex');
ylabel('$\Delta A(t)$','Interpreter','latex');
set(gca,'FontSize',32);
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

