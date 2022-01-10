% close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time = 2; lw = 2;

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end 

load('rectangle_wetting_data128_cb1e2_gamma1_0_5_reprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1; 
Nffi = para.Nffi;
pffi = reshape(hstry.pffi_track(ind,1:Nffi,1:2),Nffi,2);
kffi = getcurvature(pffi); kffi(1) = kffi(2); kffi(Nffi) = kffi(Nffi-1);
p1 = plot(pffi(1:Nffi,1),kffi(1:Nffi),':','LineWidth',lw,'color',c1);
hold on;

gap = 20;
load('rectangle_wetting_data256_cb1e_3_gamma1_0_5_reprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
Nffi = hstry.Nffis(ind);
pffi = reshape(hstry.pffi_track(ind,1:Nffi,1:2),Nffi,2);
kffi = getcurvature(pffi); kffi(1) = 2*kffi(2)-kffi(3); kffi(Nffi) = 2*kffi(Nffi-1)-kffi(Nffi-2);
p2 = plot(pffi(1+gap:Nffi-gap,1),kffi(1+gap:Nffi-gap),'-.','LineWidth',lw,'color',c1);
 
load('rectangle_dewetting_data128_cb1e2_gamma2_0_5_noprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
Nffi = para.Nffi;
pffi = reshape(hstry.pffi_track(ind,1:Nffi,1:2),Nffi,2);
kffi = getcurvature(pffi); kffi(1) = kffi(2); kffi(Nffi) = kffi(Nffi-1);
p3 = plot(pffi(1:Nffi,1),kffi(1:Nffi),'--','LineWidth',lw,'color',c2);

gap = 10;
load('rectangle_dewetting_data256_cb1e_3_gamma2_0_5_reprod_Ca_0_2.mat');
ind = floor(Time/para.dt)+1;
Nffi = hstry.Nffis(ind);
pffi = reshape(hstry.pffi_track(ind,1:Nffi,1:2),Nffi,2);
kffi = getcurvature(pffi); kffi(1) = 2*kffi(2)-kffi(3); kffi(Nffi) = 2*kffi(Nffi-1)-kffi(Nffi-2);
p4 = plot(pffi(1+gap:Nffi-gap,1),kffi(1+gap:Nffi-gap),'-','LineWidth',lw,'color',c2);

leg = legend([p1,p2,p3,p4],{'$\theta_Y=\pi/3,\ c_b=100$','$\theta_Y=\pi/3,\ c_b=0.001$','$\theta_Y=2\pi/3,\ c_b=100$','$\theta_Y=2\pi/3,\ c_b=0.001$'});
set(leg,'Interpreter','latex'); 
xlabel('$x$','Interpreter','latex');
ylabel('$\kappa$ on the fluid-fluid interface','Interpreter','latex');
set(gca,'FontSize',32);
% axis equal;
% ylim([-0.3,0.8]); 
% xlim([-1,1]);
hold off;
