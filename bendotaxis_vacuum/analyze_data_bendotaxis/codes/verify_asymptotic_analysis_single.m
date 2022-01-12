clear;
clc;
addpath(genpath('.'));
dbstop if error

load('wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1.mat');

% itr = hstry.itr; 
dt = para.dt;
itr = hstry.itr;
xr = hstry.pffir_track(1:itr,1,1);
xl = hstry.pffil_track(1:itr,1,1);

plot((0:itr-1)*dt, xr-xl, '-r');
plot((0:itr-1)*dt, xl, '-r');

% according to xr find start
for i =2:itr
%     if xl(i)>xl(i-1)
    if xr(i)>=6.14
        xr(i)
        start = i;
        itr = start + 15/dt;
        break
    end
end
% start = 100;
logt = log((1:itr-start)*dt);
% logt = log((start+1:itr)*dt);
% dxl = 1/xl(start)-1./xl(start+1:itr);
dxr = 1/xr(start)-1./xr(start+1:itr);
% logdxl = log(dxl);
logdxr = log(dxr);

plot(logt(1:itr-start),logdxr(1:itr-start),'-b','LineWidth',2.0);
% hold on;

ratior = zeros(itr,1);
for i = 3:itr-start
    ratior(i-2) = (logdxr(i-1)-logdxr(i-2))/(logt(i-1)-logt(i-2));
end
plot(logt(1:itr-start-2),ratior(1:itr-start-2),'-or');

% ratiol = zeros(itr,1);
% for i = 3:itr-start
%     ratiol(i-2) = (logdxl(i-1)-logdxl(i-2))/(logt(i-1)-logt(i-2));
% end
% plot(logt(1:itr-start-2),ratiol(1:itr-start-2),'-or');

% velxr = (xr(2:itr)-xr(1:itr-1))/dt;
% xrsquare = xr.^2;
% logvel = log(velxr);
% logxrsquare = log(xrsquare);
% 
% plot(logxrsquare(100:itr-1), logvel(100:itr-1),'-b');
% 
% ratio = zeros(itr,1);
% for i = 3:itr
%     ratio(i-2) = (logvel(i-1)-logvel(i-2))/(logxrsquare(i-1)-logxrsquare(i-2));
% end
% plot(logxrsquare(1:itr-2),ratio(1:itr-2),'-r');