close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error

% ind: 1 4 9 17
load('record_wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1_ind_17.mat'); 
load('plotvel_wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1_ind_17.mat');

disp(['Time = ',num2str((hstry.itr-1)*para.dt)]);

node = geom.node; elem = geom.elem; vel = sln.flu_vel; mem2node = geom.mem2node; mem2edge = geom.mem2edge;
N = length(node); Ne = length(elem); elem2edge = geom.elem2edge;
upbd2node = geom.upbd2node; upbd2edge = geom.upbd2edge;
pmem = sln.pmemold; Nmem = para.Nmem; By = para.By; cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2);
pffil = sln.pffilold; pffir = sln.pffirold; Nffil = para.Nffil; Nffir = para.Nffir;
pressure = sln.pressure; Nt = size(elem,1);

% [cnode,celem] = produce_coarse_mesh(geom,para,sln);
% save('plotvel_wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1_ind_9.mat','cnode','celem');
% load('plotvel_wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1_ind_9.mat');

newvel = zeros(size(cnode));

for k = 1:size(cnode,1)
    x = cnode(k,1); y = cnode(k,2); updated = false;
    for ind = 1:Ne
        tx = node(elem(ind,1:3),1); ty = node(elem(ind,1:3),2);
        if inpolygon(x,y,tx,ty) == 0 
            continue; 
        end
        ref = ([tx(2)-tx(1) tx(3)-tx(1); ty(2)-ty(1) ty(3)-ty(1)])^-1*([x-tx(1);y-ty(1)]);
        X = ref(1); Y = ref(2);
        weights = [vel(elem(ind,1:3),1:2); vel(elem2edge(ind,1:3)+N,1:2)];
        bases = [(1-X-Y)*(1-2*X-2*Y) X*(2*X-1) Y*(2*Y-1) 4*X*Y 4*Y*(1-X-Y) 4*X*(1-X-Y)];
        newvel(k,1:2) = bases*weights;
        updated = true;
        break;
    end
    if ~updated
        if abs(y-By)<1e-8
            for j = 1:length(upbd2node)-1
                if node(upbd2node(j),1)>x || node(upbd2node(j+1),1)<x
                    continue;
                end
                S = sqrt((x-node(upbd2node(j),1))^2+(y-node(upbd2node(j),2))^2)/...
                    sqrt((node(upbd2node(j+1),1)-node(upbd2node(j),1))^2+(node(upbd2node(j+1),2)-node(upbd2node(j),2))^2);
                weights = [vel(upbd2node(j:j+1),1:2); vel(upbd2edge(j)+N,1:2)];
                bases = [2*(S-1/2)*(S-1) 2*S*(S-1/2) -4*S*(S-1)];
                newvel(k,1:2) = bases*weights;
                break;
            end
        else
            for j = cl2mem1:cl2mem2-1
                if pmem(j,1)>x || pmem(j+1,1)<x
                    continue;
                end
                S = sqrt((x-pmem(j,1))^2+(y-pmem(j,2))^2)/sqrt((pmem(j+1,1)-pmem(j,1))^2+(pmem(j+1,2)-pmem(j,2))^2);
                ind = j-cl2mem1+1;
                weights = [vel(mem2node(ind:ind+1),1:2); vel(mem2edge(ind)+N,1:2)];
                bases = [2*(S-1/2)*(S-1) 2*S*(S-1/2) -4*S*(S-1)];
                newvel(k,1:2) = bases*weights;
                break;
            end  
        end
    end
end

% % trimesh(celem,cnode(:,1),cnode(:,2)-para.By,zeros(size(cnode,1),1));
% % trimesh(elem,node(:,1),node(:,2)-para.By,zeros(size(node,1),1));
% trimesh(elem,node(:,1),node(:,2)-para.By,pressure(1:N),'facecolor','interp','edgecolor','interp'); 
% hold on;
% trimesh(elem,node(:,1),para.By-node(:,2),pressure(1:N),'facecolor','interp','edgecolor','interp');
% % trimesh(celem,cnode(:,1),cnode(:,2)-para.By,zeros(size(cnode,1),1));
% % quiver(node(:,1),node(:,2)-para.By,vel(1:N,1),vel(1:N,2),'g','LineWidth',1.5);
% plot(pffil(1:Nffil,1),pffil(1:Nffil,2)-para.By,'-r','LineWidth',1.0);
% plot(pffil(1:Nffil,1),para.By-pffil(1:Nffil,2),'-r','LineWidth',1.0);
% plot(pffir(1:Nffir,1),pffir(1:Nffir,2)-para.By,'-r','LineWidth',1.0);
% plot(pffir(1:Nffir,1),para.By-pffir(1:Nffir,2),'-r','LineWidth',1.0);
% plot(pmem(1:Nmem,1),pmem(1:Nmem,2)-para.By,'-r','LineWidth',1.0);
% plot(pmem(1:Nmem,1),para.By-pmem(1:Nmem,2),'-r','LineWidth',1.0);
% plot([0 0],[-para.By, para.By],'-k','LineWidth',3.0); 
% 
% quiver(cnode(:,1),cnode(:,2)-para.By,newvel(:,1),newvel(:,2),'k','LineWidth',1.5,'AutoScaleFactor',0.5);
% quiver(cnode(:,1),para.By-cnode(:,2),newvel(:,1),-newvel(:,2),'k','LineWidth',1.5,'AutoScaleFactor',0.5);
% 
% view(2);
% axis equal;
% xlim([0, para.Bx]); ylim([-para.By-0.2, para.By+0.2]);
% colorbar('northoutside');
% set(gca,'FontSize',32);
% set(gca,'YAxisLocation','right');
% hold off;
     
gap = 10;

set(gcf,'units','Normalized','position',[0.02,0.3,0.95,0.31]);
Fig = tight_subplot(1,1,[.001 .001], [.00001 .00001], [.006 .048]);
axes(Fig(1));
trimesh(elem,node(:,1),node(:,2)-para.By,pressure(1:N),'facecolor','interp','edgecolor','interp'); 
hold on;
trimesh(elem,node(:,1),para.By-node(:,2),pressure(1:N),'facecolor','interp','edgecolor','interp');
plot(pffil(1:Nffil,1),pffil(1:Nffil,2)-para.By,'-r','LineWidth',1.0);
plot(pffil(1:Nffil,1),para.By-pffil(1:Nffil,2),'-r','LineWidth',1.0);
plot(pffir(1:Nffir,1),pffir(1:Nffir,2)-para.By,'-r','LineWidth',1.0);
plot(pffir(1:Nffir,1),para.By-pffir(1:Nffir,2),'-r','LineWidth',1.0);
% plot(pmem(cl2mem1-gap:cl2mem2+gap,1),pmem(cl2mem1-gap:cl2mem2+gap,2)-para.By,'-r','LineWidth',1.0);
% plot(pmem(cl2mem1-gap:cl2mem2+gap,1),para.By-pmem(cl2mem1-gap:cl2mem2+gap,2),'-r','LineWidth',1.0);
plot(pmem(1:Nmem,1),pmem(1:Nmem,2)-para.By,'-r','LineWidth',1.0);
plot(pmem(1:Nmem,1),para.By-pmem(1:Nmem,2),'-r','LineWidth',1.0);
plot([0 0],[-para.By, para.By],'-k','LineWidth',3.0); 

quiver(cnode(:,1),cnode(:,2)-para.By,newvel(:,1),newvel(:,2),'k','LineWidth',1.5,'AutoScaleFactor',0.6);
quiver(cnode(:,1),para.By-cnode(:,2),newvel(:,1),-newvel(:,2),'k','LineWidth',1.5,'AutoScaleFactor',0.6);

view(2);
axis equal;
% axis tight;
% xlim([pmem(cl2mem1-gap,1), pmem(cl2mem2+gap,1)]); 
xlim([0, para.Bx]);
ylim([-para.By-0.05, para.By+0.05]);
xticks([0 1 2 3 4 5 6 7 8 9 10]);
colorbar('northoutside');

% caxis([-20 -4.8]); % 1
% caxis([-17.5,-13.5]); % 4
% caxis([-17.5,-13]); % 9
caxis([-26,-13]); % 17

colorbar('northoutside','position',[0.01 0.8 0.93 0.06]);
cbarrow;

set(gca,'FontSize',32);
set(gca,'YAxisLocation','right');
hold off;

max(pressure(1:N))
min(pressure(1:N))
% elem_mid(:,1) = 1/3*(node(elem(:,1),1)+node(elem(:,2),1)+node(elem(:,3),1));
% elem_mid(:,2) = 1/3*(node(elem(:,1),2)+node(elem(:,2),2)+node(elem(:,3),2));
% % pos = [node; elem_mid];
% % plot3(pos(:,1),pos(:,2),pressure,'or'); 
% % plot3(elem_mid(:,1),elem_mid(:,2),pressure(N+1:Nt+N),'or'); 
% % plot3(node(:,1),node(:,2),pressure(1:N),'or'); 
% trimesh(elem,node(:,1),node(:,2)-para.By,pressure(1:N),'facecolor','interp','edgecolor','interp');
% % hold on;
% % view(2);
% % axis off;
% % img = frame2im(getframe(gca));
% % hold off;

