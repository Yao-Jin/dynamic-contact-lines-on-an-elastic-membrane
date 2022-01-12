clear;
clc;
addpath(genpath('.'));
dbstop if error

% Time_duration = 180; 
start_x = 6.16;

load('wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1.mat');
dt1 = para.dt; xr1 = hstry.pffir_track(1:hstry.itr,1,1);
% according to xr find start
for i =1:hstry.itr
    if xr1(i)>=start_x
        start1 = i;
        break
    end
end
end1 = start1 + 19/dt1;
% end1 = hstry.itr;
logt1 = log((1:end1-start1)*dt1);
dxr1 = 1/xr1(start1)-1./xr1(start1+1:end1); logdxr1 = log(dxr1);

load('wetting_20_10_gap_4_xl_2_gamma1_0_5_cb_5e3_Ca_0_1.mat');
xr2 = hstry.pffir_track(1:hstry.itr,1,1);
% according to xr find start
for i =1:hstry.itr
    if xr2(i)>=start_x
        start2 = i;
        break
    end
end
end2 = start2 + 35/dt1;
% end2 = hstry.itr;
logt2 = log((1:end2-start2)*dt1);
dxr2 = 1/xr2(start2)-1./xr2(start2+1:end2); logdxr2 = log(dxr2);

load('wetting_20_10_gap_4_xl_2_gamma1_0_7_cb_5e3_Ca_0_1.mat');
xr3 = hstry.pffir_track(1:hstry.itr,1,1);
% according to xr find start
for i =1:hstry.itr
    if xr3(i)>=start_x
        start3 = i;
        break
    end
end
end3 = start3 + 75/dt1;
% end3 = hstry.itr;
logt3 = log((1:end3-start3)*dt1);
dxr3 = 1/xr3(start3)-1./xr3(start3+1:end3); logdxr3 = log(dxr3);

load('wetting_20_10_gap_4_xl_2_gamma1_0_9_cb_5e3_Ca_0_1.mat');
xr4 = hstry.pffir_track(1:hstry.itr,1,1);
% according to xr find start
for i =1:hstry.itr
    if xr4(i)>=start_x
        start4 = i;
        break
    end
end
end4 = start4 + 88/dt1;
% end4 = hstry.itr;
logt4 = log((1:end4-start4)*dt1);
dxr4 = 1/xr4(start4)-1./xr4(start4+1:end4); logdxr4 = log(dxr4);

load('dewetting_20_10_gap_4_xl_2_gamma2_0_1_cb_5e3_Ca_0_1.mat');
xr5 = hstry.pffir_track(1:hstry.itr,1,1); dt2 = para.dt;
% according to xr find start
for i =2:hstry.itr
    if xr5(i)>=6.0
        start5 = i;
        break
    end
end
end5 = start5 + 60/dt2;
% end5 = hstry.itr;
logt5 = log((1:end5-start5)*dt2);
dxr5 = 1/xr5(start5)-1./xr5(start5+1:end5); logdxr5 = log(dxr5);

load('dewetting_20_10_gap_4_xl_2_gamma2_0_3_cb_5e3_Ca_0_1.mat');
xr6 = hstry.pffir_track(1:hstry.itr,1,1);
% according to xr find start
for i =2:hstry.itr
    if xr6(i)>=6.0
        start6 = i;
        break
    end
end
end6 = start6 + 88/dt2;
% end6 = end4;
logt6 = log((1:end6-start6)*dt2);
dxr6 = 1/xr6(start6)-1./xr6(start6+1:end6); logdxr6 = log(dxr6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end

% p1 = plot(logt1(2:end1-start1),logdxr1(2:end1-start1),'-','LineWidth',2.0,'color',c1);
% hold on;
% p2 = plot(logt2(2:end2-start2),logdxr2(2:end2-start2),'-.','LineWidth',2.0,'color',c1);
% p3 = plot(logt3(2:end3-start3),logdxr3(2:end3-start3),'--','LineWidth',2.0,'color',c1);
% p4 = plot(logt4(2:end4-start4),logdxr4(2:end4-start4),':','LineWidth',2.0,'color',c1);
% p5 = plot(logt5(1:end5-start5),logdxr5(1:end5-start5),'-','LineWidth',2.0,'color',c2);
% p6 = plot(logt6(1:end6-start6),logdxr6(1:end6-start6),'-.','LineWidth',2.0,'color',c2);
% plot([2 4 4 2],[-10 -10 -8 -10],'-k','LineWidth',2.0);
% 
% xlabel('$\log(t)$','Interpreter','latex');
% ylabel('$\log(\frac{1}{x_r(0)}-\frac{1}{x_r(t)})$','Interpreter','latex');
% hh = legend([p1,p2,p3,p4,p5,p6],{'$\cos\theta_Y=0.7$','$\cos\theta_Y=0.5$','$\cos\theta_Y=0.3$','$\cos\theta_Y=0.1$',...
%                 '$\cos\theta_Y=-0.9$','$\cos\theta_Y=-0.7$'});
% set(hh,'Interpreter','latex','FontName','Times New Roman','FontSize',24,'FontWeight','normal','Location','southeast');
% axis equal;
% xlim([-3,7]);
% ylim([-15,-2]);
% set(gca,'FontSize',32);
% hold off;

p1 = loglog((1:end1-start1)*dt1, dxr1(1:end1-start1),'-','LineWidth',2.0,'color',c1);
hold on;
p2 = loglog((1:end2-start2)*dt1, dxr2(1:end2-start2),'-.','LineWidth',2.0,'color',c1);
p3 = loglog((1:end3-start3)*dt1, dxr3(1:end3-start3),'--','LineWidth',2.0,'color',c1);
p4 = loglog((1:end4-start4)*dt1, dxr4(1:end4-start4),':','LineWidth',2.0,'color',c1);
p5 = loglog((1:end5-start5)*dt2, dxr5(1:end5-start5),'-','LineWidth',2.0,'color',c2);
p6 = loglog((1:end6-start6)*dt2, dxr6(1:end6-start6),'-.','LineWidth',2.0,  'color',c2);
plot([1e-1 1 1 1e-1],[1e-5 1e-5 1e-4 1e-5],'-k','LineWidth',2.0);

xlabel('$t$','Interpreter','latex');
ylabel('$\frac{1}{x_r(0)}-\frac{1}{x_r(t)}$','Interpreter','latex');
hh = legend([p1,p2,p3,p4,p5,p6],{'$\cos\theta_Y=0.7$','$\cos\theta_Y=0.5$','$\cos\theta_Y=0.3$','$\cos\theta_Y=0.1$',...
                '$\cos\theta_Y=-0.9$','$\cos\theta_Y=-0.7$'});
set(hh,'Interpreter','latex','FontName','Times New Roman','FontSize',26,'FontWeight','normal','Location','northwest');
% axis equal;
xlim([0,110]);
ylim([0 0.1]);
set(gca,'FontSize',32);
hold off;
