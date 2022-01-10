% close all;
clear;
clc;
addpath(genpath('.'));
load('initdata_adapt_32_cl_rectangle.mat');

dbstop if error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para.cb = 1e-1;
para.dt = 0.02/4^0;
para.eta1 = 1; para.mu1 = 1; para.gamma3 = 1; % fixed = 1, since dimensionless
para.Ca = 0.2; para.eta2 = 0.1; para.mu2 = 0.1; para.mu0 = 0.1; para.ls = 0.1;
para.gamma1 = 1; para.gamma2 = 0.5; % dewetting
% para.gamma1 = 0.5; para.gamma2 = 1; % wetting
para.thetaY = acos((para.gamma2-para.gamma1)/para.gamma3);

itr = 1; T = 2; maxitr = floor(T/para.dt)+1; curind = 0;

site = 'rectangle_dewetting_data32_cb1e_1_gamma2_0_5_Ca_0_2.mat';
% site = 'rectangle_wetting_data32_cb1e_1_gamma1_0_5_Ca_0_2.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------------------
% record the initial settings
for tmp = 1:1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hstry = [];
hstry.prod_mesh_times = 0;
hstry.prod_mesh_index = zeros(2,1);
hstry.memlen = zeros(maxitr,1);
hstry.droparea = zeros(maxitr,1);
hstry.Nmems = zeros(maxitr,1);
hstry.Nffis = zeros(maxitr,1);
hstry.cl2mem = zeros(maxitr,2);
hstry.pmem_track = zeros(maxitr,2*para.Nmem-1,2);
hstry.pffi_track = zeros(maxitr,2*para.Nffi-1,2);
hstry.kmem_track = zeros(maxitr,2*para.Nmem-1);
hstry.nu_track = zeros(maxitr,2*para.Nmem+1,1);
hstry.fluvel_max_norm = zeros(maxitr,1);
hstry.dynamic_angle_l = zeros(maxitr,1);
hstry.dynamic_angle_r = zeros(maxitr,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vec1 = sln.pffi(2,:)-sln.pffi(1,:);
    vec2 = sln.pmem(geom.cl2mem(1)+1,:)-sln.pmem(geom.cl2mem(1),:);
    len1 = sqrt(vec1(:,1)^2+vec1(:,2)^2);
    len2 = sqrt(vec2(:,1)^2+vec2(:,2)^2);
    thetad = acos(sum(vec1.*vec2)/(len1*len2));
    disp(['left thetad: ',num2str(thetad)]);
    hstry.dynamic_angle_l(itr) = thetad;   

    vec1 = sln.pffi(para.Nffi-1,:)-sln.pffi(para.Nffi,:);
    vec2 = sln.pmem(geom.cl2mem(2)-1,:)-sln.pmem(geom.cl2mem(2),:);
    len1 = sqrt(vec1(:,1)^2+vec1(:,2)^2);
    len2 = sqrt(vec2(:,1)^2+vec2(:,2)^2);
    thetad = acos(sum(vec1.*vec2)/(len1*len2));
    disp(['right thetad: ',num2str(thetad)]);
    hstry.dynamic_angle_r(itr) = thetad;
       
    hstry.pmem_track(itr,1:2*para.Nmem-1,1:2) = sln.pmem;
    hstry.pffi_track(itr,1:2*para.Nffi-1,1:2) = sln.pffi;
    hstry.memlen(itr) = para.initlen;
    hstry.droparea(itr) = para.areainit;
    hstry.Nmems(itr) = para.Nmem;
    hstry.Nffis(itr) = para.Nffi;
    hstry.cl2mem(itr,1:2) = geom.cl2mem;
    hstry.nodeN(itr) = size(geom.node,1);
end
% ------------------------------------------------------------------------------------
% semi-implicit scheme of the coupled system
while itr<=maxitr
    
    itr = itr + 1;
    hstry.itr = itr;
    disp('-------------------------------------------------------------------------------------------------');
    disp(['step: ',num2str(itr)]);
    
    % step 1: coupled linear system for one time step
    % update membrane
    [sln, para] = fluid_P2P01_mem_P2P2_nu_P2(geom,para,sln);
    
%     curtime = (itr-1)*para.dt; target_time = [0.02 0.04 0.05 0.06 0.1 0.2 0.5 1 2];
%     if ismember(curtime, target_time)
%         curind = curind + 1;
%         cursite = ['rectangle_wetting_data_128_ind_',num2str(curind),'.mat'];
%         save(cursite,'geom','para','sln');
%     end    

    % update fluid-fluid interface
    [sln, para] = fluid_P2P01_ffi_P2_nu_P2(geom,para,sln);

% ------------------------------------------------------------------------------------
% record tracks and disp some info
for tmp = 1:1
    vec1 = sln.pffi(2,:)-sln.pffi(1,:);
    vec2 = sln.pmem(geom.cl2mem(1)+1,:)- sln.pmem(geom.cl2mem(1),:);
    len1 = sqrt(vec1(:,1)^2+vec1(:,2)^2);
    len2 = sqrt(vec2(:,1)^2+vec2(:,2)^2);
    thetad = acos(sum(vec1.*vec2)/(len1*len2));
    disp(['left thetad: ',num2str(thetad)]);
    hstry.dynamic_angle_l(itr) = thetad;   

    vec1 = sln.pffi(para.Nffi-1,:)-sln.pffi(para.Nffi,:);
    vec2 = sln.pmem(geom.cl2mem(2)-1,:)- sln.pmem(geom.cl2mem(2),:);
    len1 = sqrt(vec1(:,1)^2+vec1(:,2)^2);
    len2 = sqrt(vec2(:,1)^2+vec2(:,2)^2);
    thetad = acos(sum(vec1.*vec2)/(len1*len2));
    disp(['right thetad: ',num2str(thetad)]);
    hstry.dynamic_angle_r(itr) = thetad;
       
    hstry.pmem_track(itr,1:2*para.Nmem-1,1:2) = sln.pmem(1:2*para.Nmem-1,1:2);
    hstry.pffi_track(itr,1:2*para.Nffi-1,1:2) = sln.pffi(1:2*para.Nffi-1,1:2);
    hstry.kmem_track(itr,1:2*para.Nmem-1) = sln.kmem(1:2*para.Nmem-1);
    hstry.nu_track(itr,1:2*para.Nmem+1) = sln.nu(1:2*para.Nmem+1);
    hstry.memlen(itr) = para.curlen;
    hstry.droparea(itr) = para.areanow;
    hstry.Nmems(itr) = para.Nmem;
    hstry.Nffis(itr) = para.Nffi;
    hstry.cl2mem(itr,1:2) = geom.cl2mem;
    hstry.fluvel_max_norm(itr) = max(sqrt(sum(sln.flu_vel.^2,2)));
end

% -----------------------------------------------------------------------------------
    
    % step 2: redistribute markers on ffi and mem
    [sln, para, reproduce_mesh1] = adaptffi_even_P2(sln,para);
    [geom, para, sln, reproduce_mesh2] = adaptmem_even_P2(geom,para,sln);
    
    % step 3: update the triangular mesh by moving mesh method or just produce a new mesh
%     geom = Movingmesh_lrperiodic(geom,para,sln);
    if (reproduce_mesh1==0) && (reproduce_mesh2==0)
        geom = Movingmesh_lrperiodic(geom,para,sln);
    else
        hstry.prod_mesh_times = hstry.prod_mesh_times+1;
        disp(['reproduce mesh times: ',num2str(hstry.prod_mesh_times)]);
        geom = adaptLRbd(geom,para,sln);
        geom = produce_newmesh(geom,para,sln);
     
        hstry.nodeN(hstry.prod_mesh_times+1) = size(geom.node,1);
        hstry.prod_mesh_index(hstry.prod_mesh_times) = itr;
    end

    if mod(itr,50)==0
       save(site,'geom','para','sln','hstry');
    end
    
end
save(site,'geom','para','sln','hstry');