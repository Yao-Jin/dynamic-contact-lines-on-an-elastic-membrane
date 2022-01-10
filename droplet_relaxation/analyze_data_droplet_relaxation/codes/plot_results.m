clear;
% clc;
addpath(genpath('.'));
dbstop if error

load('rectangle_dewetting_data256_cb1e_3_gamma2_0_5_reprod_Ca_0_2.mat');

node = geom.node; elem = geom.elem;
pos = [node; 0.5*(node(geom.edge(:,1),1:2)+node(geom.edge(:,2),1:2))];
trimesh(elem,node(:,1),node(:,2),zeros(size(node,1),1));
hold on;
%         quiver(node(:,1),node(:,2),flu_vel(1:N,1),flu_vel(1:N,2),'r');
quiver(pos(:,1),pos(:,2),sln.flu_vel(:,1),sln.flu_vel(:,2),'r');
view(2);
axis equal;
hold off;

% [sln] = interploate_fluvel(geom,para,sln);

itr = hstry.itr;
plot(para.dt*(0:itr-1),hstry.dynamic_angle_l(1:itr),'-or');
hold on;
plot([0 para.dt*(itr-1)],[para.thetaY para.thetaY],'-k');
hold off;

plot(para.dt*(0:itr-1),hstry.pffi_track(1:itr,1,1),'-or');

plot(para.dt*(0:itr-1),hstry.pmem_track(1:itr,1,2),'-or');


plot(para.dt*(0:itr-1),hstry.fluvel_max_norm(1:itr,1,1),'-or');

% plot(para.dt*(0:itr-1),hstry.maxdy(1:itr,1,1),'-or');

plot(para.dt*(0:itr-1),hstry.droparea(1:itr),'-or');
hold on;
plot([0 para.dt*(itr-1)],[para.areainit para.areainit],'-k');
hold off;

plot(para.dt*(0:itr-1),hstry.memlen(1:itr),'-or');


Nmem = para.Nmem; cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2);
x = sln.pmem(1:Nmem,1); y = sln.pmem(1:Nmem,2);
plot(x(1:cl2mem1),sln.bargamma(1:cl2mem1),'-ob');
hold on;
plot(x(cl2mem1:cl2mem2),[sln.bargamma(Nmem+1);sln.bargamma(cl2mem1+1:cl2mem2-1);sln.bargamma(Nmem+2)],'-sr');
plot(x(cl2mem2:Nmem),sln.bargamma(cl2mem2:Nmem),'-ob');  
hold off; 

disp(['left bargamma jump: ', num2str(sln.bargamma(Nmem+1)-sln.bargamma(cl2mem1))]);    
disp(['right bargamma jump: ', num2str(sln.bargamma(Nmem+2)-sln.bargamma(cl2mem2))]); 

kappa = getcurvature(sln.pmem(1:Nmem,:));
plot(x(1:Nmem),kappa,'-*r');

% kmem = sln.kmem(1:Nmem);
% plot(x(1:Nmem),kmem,'-*r');
% hold on;
% plot(sln.pmem(Nmem+1:end,1),sln.kmem(Nmem+1:end),'-ob');
% hold off;

