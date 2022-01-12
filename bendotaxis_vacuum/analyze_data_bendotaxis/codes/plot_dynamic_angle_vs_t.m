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

load('testconv_10_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');
ind = floor(Time/para.dt)+1;
nus = hstry.nu_track; Nm = para.Nm;
dnusL = nus(:,1)-nus(:,Nm+1);
dnusR = nus(:,Nm);
angleL = acos(para.gamma2-para.gamma1+dnusL);
angleR = acos(para.gamma2-para.gamma1+dnusR);
plot(para.dt*(0:ind-1),angleR(1:ind),':k','LineWidth',2.0);
% angleL(ind)-para.thetaY
hold on;

load('testconv_20_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');
ind = floor(Time/para.dt)+1;
nus = hstry.nu_track; Nm = para.Nm;
dnusL = nus(:,1)-nus(:,Nm+1);
dnusR = nus(:,Nm);
angleL = acos(para.gamma2-para.gamma1+dnusL);
angleR = acos(para.gamma2-para.gamma1+dnusR);
plot(para.dt*(0:ind-1),angleR(1:ind),'--','LineWidth',2.0,'color',c2);
% angleL(ind)-para.thetaY

load('testconv_40_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');
ind = floor(Time/para.dt)+1;
nus = hstry.nu_track; Nm = para.Nm;
dnusL = nus(:,1)-nus(:,Nm+1);
dnusR = nus(:,Nm);
angleL = acos(para.gamma2-para.gamma1+dnusL);
angleR = acos(para.gamma2-para.gamma1+dnusR);
plot(para.dt*(0:ind-1),angleR(1:ind),'-.','LineWidth',2.0,'color',c3);
% angleL(ind)-para.thetaY

load('testconv_80_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');
ind = floor(Time/para.dt)+1;
nus = hstry.nu_track; Nm = para.Nm;
dnusL = nus(:,1)-nus(:,Nm+1);
dnusR = nus(:,Nm);
angleL = acos(para.gamma2-para.gamma1+dnusL);
angleR = acos(para.gamma2-para.gamma1+dnusR);
plot(para.dt*(0:ind-1),angleR(1:ind),'-','LineWidth',2.0,'color',c1);
% angleL(ind)-para.thetaY

% plot([0 Time],[para.thetaY para.thetaY],'-','LineWidth',2.0,'color',c4);
leg = legend('$\mathcal{R}=10$','$\mathcal{R}=20$', '$\mathcal{R}=40$', '$\mathcal{R}=80$');
set(leg,'Interpreter','latex'); 
xlabel('$t$','Interpreter','latex');
ylabel('$\theta_d$','Interpreter','latex');
set(gca,'FontSize',32);
% ylim([1.55,2.1]);
% axis tight;
% axis equal;
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

