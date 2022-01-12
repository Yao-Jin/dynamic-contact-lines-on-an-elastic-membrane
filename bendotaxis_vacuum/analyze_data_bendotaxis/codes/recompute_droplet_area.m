close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error

site = 'testconv_80_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat';
load(site);
Nffil = para.Nffil; Nffir = para.Nffir; Nm = para.Nm;
for i = 1:hstry.itr
   pffil = reshape(hstry.pffil_track(i,1:Nffil,1:2),Nffil,2); 
   pffir = reshape(hstry.pffir_track(i,1:Nffir,1:2),Nffir,2); 
   cl2mem1 = hstry.cl2mem(i,1); cl2mem2 = hstry.cl2mem(i,2);
   pmem = reshape(hstry.pmem_track(i,cl2mem1:cl2mem2,1:2),Nm,2);
   dropletbd = [pffil;pffir(end:-1:1,1:2);pmem(Nm-1:-1:2,1:2)];
   hstry.droparea(i) = polyarea(dropletbd(:,1),dropletbd(:,2));
end
save(site,'geom','para','sln','hstry');