% close all;
clear;
clc;
addpath(genpath('.'));

load('init_vacuum_predewetting_20_10_gap_4_xl_2_gamma2_0_3.mat');
% load('init_vacuum_prewetting_20_10_gap_4_xl_2_gamma1_0_3.mat');
% load('init_vacuum_testconv_prewetting_20_3_0_5_gap_1_xl_1_gamma1_0_5.mat');

dbstop if error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para.cb = 1e2; 
para.dt = 0.02/4^0;
para.Ca = 0.1; para.mu0 = 0.1; para.ls = 0.1; para.eta1 = 1; para.mu1 = 1;
% gamma1, gamma2, gamma3, thetaY are pre-determined

itr = 1; T = 3; maxitr = floor(T/para.dt)+1;

% site = 'wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1.mat';
% site = 'dewetting_20_10_gap_4_xl_2_gamma2_0_3_cb_3e3_Ca_0_1.mat';
site = 'testconv_10_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat';

index = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------------------
% record the initial settings
for tmp = 1:1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hstry = [];
hstry.prod_mesh_times = 0;
hstry.memlen = zeros(maxitr,1);
hstry.droparea = zeros(maxitr,1);
hstry.Nmems = zeros(maxitr,1);
hstry.Nffils = zeros(maxitr,1);
hstry.Nffirs = zeros(maxitr,1);
hstry.Nms = zeros(maxitr,1);
hstry.cl2mem = zeros(maxitr,2);
hstry.pmem_track = zeros(maxitr,2*para.Nmem-1,2);
hstry.kmem_track = zeros(maxitr,2*para.Nmem-1,1);
hstry.pffil_track = zeros(maxitr,2*para.Nffil-1,2);
hstry.pffir_track = zeros(maxitr,2*para.Nffir-1,2);
hstry.nu_track = zeros(maxitr,para.Nmem,1);
hstry.dynamic_angle_l = zeros(maxitr,1);
hstry.dynamic_angle_r = zeros(maxitr,1);
hstry.fluvel_max_norm = zeros(maxitr,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vec1 = sln.pffil(2,:)-sln.pffil(1,:);
    vec2 = sln.pmem(geom.cl2mem(1)+1,:)-sln.pmem(geom.cl2mem(1),:);
    len1 = sqrt(vec1(:,1)^2+vec1(:,2)^2);
    len2 = sqrt(vec2(:,1)^2+vec2(:,2)^2);
    thetad = acos(sum(vec1.*vec2)/(len1*len2));
    disp(['left thetad: ',num2str(thetad)]);
    hstry.dynamic_angle_l(itr) = thetad;   

    vec1 = sln.pffir(2,:)-sln.pffir(1,:);
    vec2 = sln.pmem(geom.cl2mem(2)-1,:)-sln.pmem(geom.cl2mem(2),:);
    len1 = sqrt(vec1(:,1)^2+vec1(:,2)^2);
    len2 = sqrt(vec2(:,1)^2+vec2(:,2)^2);
    thetad = acos(sum(vec1.*vec2)/(len1*len2));
    disp(['right thetad: ',num2str(thetad)]);
    hstry.dynamic_angle_r(itr) = thetad;
       
    hstry.pmem_track(itr,1:2*para.Nmem-1,1:2) = sln.pmem(1:2*para.Nmem-1,1:2);
    hstry.pffil_track(itr,1:2*para.Nffil-1,1:2) = sln.pffil(1:2*para.Nffil-1,1:2);
    hstry.pffir_track(itr,1:2*para.Nffir-1,1:2) = sln.pffir(1:2*para.Nffir-1,1:2);
    hstry.memlen(itr) = para.initlen;
    hstry.droparea(itr) = para.areainit;
    hstry.Nmems(itr) = para.Nmem;
    hstry.Nffils(itr) = para.Nffil;
    hstry.Nffirs(itr) = para.Nffir;
    hstry.Nms(itr) = para.Nm;
    hstry.cl2mem(itr,1:2) = geom.cl2mem;
end
% ------------------------------------------------------------------------------------
% semi-implicit scheme of the coupled system
while itr<=maxitr
    
    itr = itr + 1;
    hstry.itr = itr;
    disp('-------------------------------------------------------------------------------------------------');
    disp(['step: ',num2str(itr)]);
    
    % step 1: coupled linear system for one time step
    % variables: u, p, y, kmem, nu
    % FEM orders: P2, P0, P2, P2, P1

    [sln, para] = fluid_P2P01_mem_P2P2_nu_P2(geom,para,sln);
    
%     curtime = (itr-1)*para.dt; target_time = [0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 10 12 15]; 
%     if ismember(curtime, target_time)
%         index = index + 1;
%         cursite = ['record_wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1_ind_',num2str(index),'.mat'];
%         cursite = ['record_wetting_20_10_gap_4_xl_2_gamma1_0_3_cb_5e3_Ca_0_1_t_0_4.mat'];
%         save(cursite,'geom','para','sln','itr');
%     end    
    
    % update ffi and contact lines   
    [sln, para] = fluid_P2P01_ffi_P2_nu_P2(geom, para, sln);  

% ------------------------------------------------------------------------------------
% record tracks and disp some info
for tmp = 1:1
    vec1 = sln.pffil(2,:)-sln.pffil(1,:);
    vec2 = sln.pmem(geom.cl2mem(1)+1,:)- sln.pmem(geom.cl2mem(1),:);
    len1 = sqrt(vec1(:,1)^2+vec1(:,2)^2);
    len2 = sqrt(vec2(:,1)^2+vec2(:,2)^2);
    thetad = acos(sum(vec1.*vec2)/(len1*len2));
    disp(['left thetad: ',num2str(thetad)]);
    hstry.dynamic_angle_l(itr) = thetad;   

    vec1 = sln.pffir(2,:)-sln.pffir(1,:);
    vec2 = sln.pmem(geom.cl2mem(2)-1,:)- sln.pmem(geom.cl2mem(2),:);
    len1 = sqrt(vec1(:,1)^2+vec1(:,2)^2);
    len2 = sqrt(vec2(:,1)^2+vec2(:,2)^2);
    thetad = acos(sum(vec1.*vec2)/(len1*len2));
    disp(['right thetad: ',num2str(thetad)]);
    hstry.dynamic_angle_r(itr) = thetad;
       
    hstry.pmem_track(itr,1:2*para.Nmem-1,1:2) = sln.pmem(1:2*para.Nmem-1,1:2);
    hstry.kmem_track(itr,1:2*para.Nmem-1) = sln.kmem(1:2*para.Nmem-1);    
    hstry.pffil_track(itr,1:2*para.Nffil-1,1:2) = sln.pffil(1:2*para.Nffil-1,1:2);
    hstry.pffir_track(itr,1:2*para.Nffir-1,1:2) = sln.pffir(1:2*para.Nffir-1,1:2);
    hstry.nu_track(itr,1:2*para.Nm) = sln.nu;
    hstry.memlen(itr) = para.curlen;
    hstry.droparea(itr) = para.areanow;
    hstry.Nmems(itr) = para.Nmem;
    hstry.Nffils(itr) = para.Nffil;
    hstry.Nffirs(itr) = para.Nffir;
    hstry.Nms(itr) = para.Nm;
    hstry.cl2mem(itr,1:2) = geom.cl2mem;
    hstry.fluvel_max_norm(itr) = max(sqrt(sum(sln.flu_vel.^2,2)));
end
% -----------------------------------------------------------------------------------
    
    % step 2: redistribute markers on mem and ffi
    [geom, para, sln] = adaptmem_even_keep(geom,para,sln);
    [sln, para] = adaptffi_even_P2(sln,para);
    
    % step 3: update the triangular mesh by moving mesh method or just produce a new mesh
    geom = Movingmesh_lrfree(geom,para,sln);

    if mod(itr,50)==0
       save(site,'geom','para','sln','hstry');
    end
    
    if abs(sln.pmem(end,2)-para.By)<0.05
        save(site,'geom','para','sln','hstry');
        break;
    end
end
save(site,'geom','para','sln','hstry');