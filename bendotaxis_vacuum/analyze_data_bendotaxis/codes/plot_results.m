clear;
% clc;
addpath(genpath('.'));
dbstop if error

load('wetting_20_10_gap_4_xl_2_gamma1_0_5_cb_5e3_Ca_0_1.mat');

N = size(geom.node,1);
trimesh(geom.elem, geom.node(:,1),geom.node(:,2),zeros(size(geom.node,1),1));
hold on;
plot(geom.node(geom.bdnode(:)==1,1),geom.node(geom.bdnode(:)==1,2),'bo');
view(2);
quiver(geom.node(:,1),geom.node(:,2),sln.flu_vel(1:N,1),sln.flu_vel(1:N,2),'r');
mid = 0.5*[geom.node(geom.edge(:,1),1)+geom.node(geom.edge(:,2),1) geom.node(geom.edge(:,1),2)+geom.node(geom.edge(:,2),2)];
quiver(mid(:,1),mid(:,2),sln.flu_vel(N+1:end,1),sln.flu_vel(N+1:end,2),'g');
axis equal;
hold off;

itr = hstry.itr;
plot(para.dt*(0:itr-1),hstry.dynamic_angle_l(1:itr),'-or');
hold on;
plot([0 para.dt*(itr-1)],[para.thetaY para.thetaY],'-k');
hold off;

plot(para.dt*(0:itr-1),hstry.fluvel_max_norm(1:itr,1,1),'-or');

plot(para.dt*(0:itr-1),hstry.droparea(1:itr),'-or');
hold on;
plot([0 para.dt*(itr-1)],[para.areainit para.areainit],'-k');
hold off;

plot(para.dt*(0:itr-1),hstry.memlen(1:itr),'-or');


Nmem = para.Nmem; cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2);
x = sln.pmem(1:Nmem,1); y = sln.pmem(1:Nmem,2); Nm = para.Nm;
plot(x(cl2mem1:cl2mem2),sln.nu(1:Nm),'-sr');
hold on;
plot([x(1) x(cl2mem1)],[sln.nu(Nm+1) sln.nu(Nm+1)],'-ob');
plot([x(cl2mem2) x(Nmem)],[0 0],'-ob');
hold off; 

disp(['left nu jump: ', num2str(sln.nu(1)-sln.nu(Nm+1))]);    
disp(['right nu jump: ', num2str(sln.nu(Nm))]); 

kmem = sln.kmem(1:Nmem);
plot(x(1:Nmem),kmem,'-*r');
hold on;
plot(sln.pmem(Nmem+1:end,1),sln.kmem(Nmem+1:end),'-ob');
hold off;

