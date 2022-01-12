% close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error

load('testconv_dewetting_64_4_gap_2_xl_1_gamma2_0_3_cb_1e1.mat');
% load('wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3.mat');

Nffil = para.Nffil; Nffir = para.Nffir; pmems = hstry.pmem_track; 
pffils = hstry.pffil_track; pffirs = hstry.pffir_track;
Bx = para.Bx; By = para.By;

i = 500; gap = 30;
Nmem = hstry.Nmems(i);
plot(pmems(i,1:Nmem,1),pmems(i,1:Nmem,2),'-b');
hold on;
plot([0 Bx], [By By],'--k');
plot([0 0], [0 By],'--k');
plot([pmems(i,Nmem-gap,1) pmems(i,Nmem-gap,1)], [pmems(i,Nmem-gap,2) By],'--k');
plot(pffils(i,1:Nffil,1),pffils(i,1:Nffil,2),'-r');
plot(pffirs(i,1:Nffir,1),pffirs(i,1:Nffir,2),'-r');

plot(pmems(i,1:Nmem,1),2*By-pmems(i,1:Nmem,2),'-b');
plot(pffils(i,1:Nffil,1),2*By-pffils(i,1:Nffil,2),'-r');
plot(pffirs(i,1:Nffir,1),2*By-pffirs(i,1:Nffir,2),'-r');
plot([0 0], [By 2*By],'--k');
plot([pmems(i,Nmem-gap,1) pmems(i,Nmem-gap,1)], [By 2*By-pmems(i,Nmem-gap,2)],'--k');


axis equal; axis off;


xmem = reshape(pmems(i,1:Nmem,1),Nmem,1);
ymem = reshape(pmems(i,1:Nmem,2),Nmem,1);
ind1 = floor(Nmem/2);
Dy = (ymem(ind1+1)-ymem(ind1-1))/(xmem(ind1+1)-xmem(ind1-1));
loclen = sqrt(1+Dy^2);
nmem1 = [Dy -1]/loclen;
tmem1 = [1 Dy]/loclen;

ind2m = floor(Nmem/8);
Dy = (ymem(ind2m+1)-ymem(ind2m-1))/(xmem(ind2m+1)-xmem(ind2m-1));
loclen = sqrt(1+Dy^2);
nmem2 = [Dy -1]/loclen;
tmem2 = [1 Dy]/loclen;

ind3m = floor(5*Nmem/6);
Dy = (ymem(ind3m+1)-ymem(ind3m-1))/(xmem(ind3m+1)-xmem(ind3m-1));
loclen = sqrt(1+Dy^2);
nmem3 = [Dy -1]/loclen;
tmem3 = [1 Dy]/loclen;


xffil = reshape(pffils(i,1:Nffil,1),Nffil,1);
yffil = reshape(pffils(i,1:Nffil,2),Nffil,1);
ind2 = floor(Nffil/2);
Dy = (yffil(ind2+1)-yffil(ind2-1))/(xffil(ind2+1)-xffil(ind2-1));
loclen = sqrt(1+Dy^2);
nffil = [-Dy 1]/loclen;
tffil = [1 Dy]/loclen;

xffir = reshape(pffirs(i,1:Nffir,1),Nffir,1);
yffir = reshape(pffirs(i,1:Nffir,2),Nffir,1);
ind3 = floor(Nffir/2);
Dy = (yffir(ind3+1)-yffir(ind3-1))/(xffir(ind3+1)-xffir(ind3-1));
loclen = sqrt(1+Dy^2);
nffir = [-Dy 1]/loclen;
tffir = [-1 -Dy]/loclen;

m_w = [xmem(Nmem-gap)-xmem(Nmem-gap-1); ymem(Nmem-gap)-ymem(Nmem-gap-1)]./...
        (sqrt((xmem(Nmem-gap)-xmem(Nmem-gap-1))^2+(ymem(Nmem-gap)-ymem(Nmem-gap-1))^2));

% pos = [xmem(ind1) xmem(ind1) xmem(ind2m) xmem(ind2m) xmem(ind3m) xmem(ind3m) xffil(ind2) xffil(ind2) xffir(ind3) xffir(ind3) 0    Bx/3 Bx    Bx/3 xmem(end); 
%        ymem(ind1) ymem(ind1) ymem(ind2m) ymem(ind2m) ymem(ind3m) ymem(ind3m) yffil(ind2) yffil(ind2) yffir(ind3) yffir(ind3) By/2 By   By/4  By   ymem(end)]';
% vec = 0.15*[nmem1; tmem1; nmem2; tmem2; nmem3; tmem3; -nffil; -tffil; -nffir; -tffir; -1 0; 0 1; 1 0; 1 0; m_w(1) m_w(2)];

pos = [xmem(ind1) xmem(ind1) xffil(ind2) xffil(ind2) xffir(ind3) xffir(ind3) 0    Bx/3 xmem(Nmem-gap) Bx/3 xmem(Nmem-gap); 
       ymem(ind1) ymem(ind1) yffil(ind2) yffil(ind2) yffir(ind3) yffir(ind3) By/2 By   By/4  By   ymem(Nmem-gap)]';
vec = 0.15*[nmem1; tmem1; -nffil; -tffil; -nffir; -tffir; -1 0; 0 1; 1 0; 1 0; m_w(1) m_w(2)];

quiver(pos(:,1),pos(:,2),vec(:,1),vec(:,2),'-k','MaxHeadSize',1,'LineWidth',1,'AutoScale','off');

plot([xmem(Nmem-gap) xmem(Nmem-gap)],[ymem(Nmem-gap) ymem(Nmem-gap)-0.4],'k--');
plot([Bx Bx+0.5],[By By],'k--');
plot([0 -0.5],[By By],'k--');
plot([0 0],[0 -0.5],'k--');
plot([0 -0.5],[0 0],'k--');

% annotation('textbox',[.1 .1 .5 .5],'String','$\mathbf{n}$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% annotation('textbox',[.15 .1 .5 .5],'String','$\mathbf{n}_w$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% annotation('textbox',[.1 .1 .6 .6],'String','\boldmath{$\tau$}','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% annotation('textbox',[.1 .5 .5 .5],'String','$\Sigma_L$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% annotation('textbox',[.9 .1 .5 .5],'String','$\Sigma_R$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% annotation('textbox',[.5 .5 .5 .5],'String','$\Sigma_U$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% annotation('textbox',[.6 .5 .5 .5],'String','$\Omega_1$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% annotation('textbox',[.8 .5 .5 .5],'String','$\Omega_2$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% annotation('textbox',[.9 .5 .5 .5],'String','$\Omega_2$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% annotation('textbox',[.7 .5 .5 .5],'String','$\Sigma_1$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% annotation('textbox',[.2 .5 .5 .5],'String','$\Sigma_2$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% annotation('textbox',[.3 .5 .5 .5],'String','$\Sigma_2$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% annotation('textbox',[.4 .5 .5 .5],'String','$\Sigma_3$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% annotation('textbox',[.45 .5 .5 .5],'String','$\Sigma_4$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% 
% annotation('textbox',[.66 .5 .5 .5],'String','$\Lambda_l$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');
% annotation('textbox',[.55 .5 .5 .5],'String','$\Lambda_r$','FontSize',20,'FitBoxToText','on','LineStyle','none','interpreter','latex');


hold off;
