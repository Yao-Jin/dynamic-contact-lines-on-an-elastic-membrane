close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error

% load('testconv_dewetting_16_4_gap_2_xl_1_gamma2_0_3_cb_1e2.mat');
% load('testconv_16_cb_1e2.mat');
% load('testconv4_10_cb_1e2.mat');
load('testconv_10_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');

Nmem1s = hstry.Nmems; Nffil1s = hstry.Nffils; Nffir1s = hstry.Nffirs;
pmem1s = hstry.pmem_track; pffil1s = hstry.pffil_track; pffir1s = hstry.pffir_track;

% load('testconv_dewetting_32_4_gap_2_xl_1_gamma2_0_3_cb_1e2.mat');
% load('testconv_32_cb_1e2.mat');
% load('testconv4_20_cb_1e2.mat');
load('testconv_20_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');

Nmem2s = hstry.Nmems; Nffil2s = hstry.Nffils; Nffir2s = hstry.Nffirs;
pmem2s = hstry.pmem_track; pffil2s = hstry.pffil_track; pffir2s = hstry.pffir_track;

% load('testconv_dewetting_64_4_gap_2_xl_1_gamma2_0_3_cb_1e2.mat');
% load('testconv_64_cb_1e2.mat');
% load('testconv4_40_cb_1e2.mat');
load('testconv_40_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');

Nmem3s = hstry.Nmems; Nffil3s = hstry.Nffils; Nffir3s = hstry.Nffirs;
pmem3s = hstry.pmem_track; pffil3s = hstry.pffil_track; pffir3s = hstry.pffir_track;

load('testconv_80_3_0_5_gap_1_xl_1_gamma1_0_5_cb_1e2_Ca_0_1.mat');

Nmem4s = hstry.Nmems; Nffil4s = hstry.Nffils; Nffir4s = hstry.Nffirs;
pmem4s = hstry.pmem_track; pffil4s = hstry.pffil_track; pffir4s = hstry.pffir_track;

aa = 150; original_dt = 0.02;
order1 = zeros(aa,1); order2 = zeros(aa,1);
for itr = 2:1:aa+1

ind = itr;
Nmem1 = Nmem1s(ind); Nffil1 = Nffil1s(ind); Nffir1 = Nffir1s(ind);
pmem1 = reshape(pmem1s(ind,1:Nmem1,1:2),Nmem1,2);
pffil1 = reshape(pffil1s(ind,1:Nffil1,1:2),Nffil1,2);
pffir1 = reshape(pffir1s(ind,1:Nffir1,1:2),Nffir1,2);

ind = (itr-1)*4+1;
Nmem2 = Nmem2s(ind); Nffil2 = Nffil2s(ind); Nffir2 = Nffir2s(ind);
pmem2 = reshape(pmem2s(ind,1:Nmem2,1:2),Nmem2,2);
pffil2 = reshape(pffil2s(ind,1:Nffil2,1:2),Nffil2,2);
pffir2 = reshape(pffir2s(ind,1:Nffir2,1:2),Nffir2,2);

ind = (itr-1)*16+1;
Nmem3 = Nmem3s(ind); Nffil3 = Nffil3s(ind); Nffir3 = Nffir3s(ind);
pmem3 = reshape(pmem3s(ind,1:Nmem3,1:2),Nmem3,2);
pffil3 = reshape(pffil3s(ind,1:Nffil3,1:2),Nffil3,2);
pffir3 = reshape(pffir3s(ind,1:Nffir3,1:2),Nffir3,2);

ind = (itr-1)*64+1;
Nmem4 = Nmem4s(ind); Nffil4 = Nffil4s(ind); Nffir4 = Nffir4s(ind);
pmem4 = reshape(pmem4s(ind,1:Nmem4,1:2),Nmem4,2);
pffil4 = reshape(pffil4s(ind,1:Nffil4,1:2),Nffil4,2);
pffir4 = reshape(pffir4s(ind,1:Nffir4,1:2),Nffir4,2);

error1 = 0;
for i = 1:Nmem1
    dist = abs(p_poly_dist(pmem1(i,:),pmem2));
    if dist>error1
        error1 = dist;
    end
end

for i = 1:Nffil1
    dist = abs(p_poly_dist(pffil1(i,:),pffil2));
    error1 = max(error1, dist);
end

for i = 1:Nffir1
    dist = abs(p_poly_dist(pffir1(i,:),pffir2));
    error1 = max(error1, dist);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error2 = 0;
for i = 1:Nmem2
    dist = abs(p_poly_dist(pmem2(i,:),pmem3));
    if dist>error2
        error2 = dist;
    end    
end

for i = 1:Nffil2
    dist = abs(p_poly_dist(pffil2(i,:),pffil3));
    error2 = max(error2, dist);
end

for i = 1:Nffir2
    dist = abs(p_poly_dist(pffir2(i,:),pffir3));
    error2 = max(error2, dist);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error3 = 0;
for i = 1:Nmem3
    dist = abs(p_poly_dist(pmem3(i,:),pmem4));
    if dist>error3
        error3 = dist;
    end    
end

for i = 1:Nffil3
    dist = abs(p_poly_dist(pffil3(i,:),pffil4));
    error3 = max(error3, dist);
end

for i = 1:Nffir3
    dist = abs(p_poly_dist(pffir3(i,:),pffir4));
    error3 = max(error3, dist);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

order1(itr-1) = log(error1/error2)/log(2);
order2(itr-1) = log(error2/error3)/log(2);

if ismember((itr-1)*original_dt,[0.4,1,2,3])
    error1
    error2
    error3
    order1(itr-1)
    order2(itr-1)
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot((1:aa)*original_dt,order1,'-or');
hold on;
plot((1:aa)*original_dt,order2,'-sb');
hold off;