clear;
clc;
addpath(genpath('.'));
dbstop if error

Time_duration = 180; start_x = 6.16;

load('wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1.mat');
dt1 = para.dt; xr1 = hstry.pffir_track(1:hstry.itr,1,1); xl1 = hstry.pffil_track(1:hstry.itr,1,1);
% according to xr find start
for i =1:hstry.itr
    if xr1(i)>=start_x
        start1 = i;
        break
    end
end
end1 = start1 + 19/dt1;
% end1 = hstry.itr;
dx1 = xr1(start1:end1)-xl1(start1:end1);

load('wetting_20_10_gap_4_xl_2_gamma1_0_5_cb_5e3_Ca_0_1.mat');
xr2 = hstry.pffir_track(1:hstry.itr,1,1); xl2 = hstry.pffil_track(1:hstry.itr,1,1);
% according to xr find start
for i =1:hstry.itr
    if xr2(i)>=start_x
        start2 = i;
        break
    end
end
end2 = start2 + 35/dt1;
% end2 = hstry.itr;
dx2 = xr2(start2:end2)-xl2(start2:end2);

load('wetting_20_10_gap_4_xl_2_gamma1_0_7_cb_5e3_Ca_0_1.mat');
xr3 = hstry.pffir_track(1:hstry.itr,1,1); xl3 = hstry.pffil_track(1:hstry.itr,1,1);
% according to xr find start
for i =1:hstry.itr
    if xr3(i)>=start_x
        start3 = i;
        break
    end
end
end3 = start3 + 75/dt1;
% end3 = hstry.itr;
dx3 = xr3(start3:end3)-xl3(start3:end3);

load('wetting_20_10_gap_4_xl_2_gamma1_0_9_cb_5e3_Ca_0_1.mat');
xr4 = hstry.pffir_track(1:hstry.itr,1,1); xl4 = hstry.pffil_track(1:hstry.itr,1,1);
% according to xr find start
for i =1:hstry.itr
    if xr4(i)>=start_x
        start4 = i;
        break
    end
end
end4 = start4 + 88/dt1;
% end4 = hstry.itr;
dx4 = xr4(start4:end4)-xl4(start4:end4);

load('dewetting_20_10_gap_4_xl_2_gamma2_0_1_cb_5e3_Ca_0_1.mat'); dt2 = para.dt;
xr5 = hstry.pffir_track(1:hstry.itr,1,1); xl5 = hstry.pffil_track(1:hstry.itr,1,1);
% according to xr find start
for i =2:hstry.itr
    if xr5(i)>=6
        start5 = i;
        break
    end
end
end5 = start5 + 60/dt2;
% end5 = hstry.itr;
dx5 = xr5(start5:end5)-xl5(start5:end5);

load('dewetting_20_10_gap_4_xl_2_gamma2_0_3_cb_5e3_Ca_0_1.mat');
xr6 = hstry.pffir_track(1:hstry.itr,1,1); xl6 = hstry.pffil_track(1:hstry.itr,1,1);
% according to xr find start
for i =2:hstry.itr
    if xr6(i)>=6
        start6 = i;
        break
    end
end
end6 = start6 + 88/dt2;
% end6 = hstry.itr;
dx6 = xr6(start6:end6)-xl6(start6:end6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end

p1 = plot((0:end1-start1)*dt1,dx1-dx1(1),'-','LineWidth',2.0,'color',c1);
hold on;
p2 = plot((0:end2-start2)*dt1,dx2-dx2(1),'-.','LineWidth',2.0,'color',c1);
p3 = plot((0:end3-start3)*dt1,dx3-dx3(1),'--','LineWidth',2.0,'color',c1);
p4 = plot((0:end4-start4)*dt1,dx4-dx4(1),':','LineWidth',2.0,'color',c1);
p5 = plot((0:end5-start5)*dt2,dx5-dx5(1),'-','LineWidth',2.0,'color',c2);
p6 = plot((0:end6-start6)*dt2,dx6-dx6(1),'-.','LineWidth',2.0,'color',c2);

hh = legend([p1,p2,p3,p4,p5,p6],{'$\cos\theta_Y=0.7$','$\cos\theta_Y=0.5$','$\cos\theta_Y=0.3$','$\cos\theta_Y=0.1$',...
                '$\cos\theta_Y=-0.9$','$\cos\theta_Y=-0.7$'});
set(hh,'Interpreter','latex','FontName','Times New Roman','FontSize',24,'FontWeight','normal','Location','northeast');

xlabel('$t$','Interpreter','latex');
ylabel('$\Delta L(t)$','Interpreter','latex');
set(gca,'FontSize',32);
% xlim([0,180]);
% ylim([-15,-3]);
hold off;
