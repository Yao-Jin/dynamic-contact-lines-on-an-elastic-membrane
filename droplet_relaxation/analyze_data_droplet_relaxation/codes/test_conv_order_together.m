close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error

load('rectangle_dewetting_data32_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');

Nmem1s = hstry.Nmems; Nffi1 = para.Nffi;
pmem1s = hstry.pmem_track;
pffi1s = hstry.pffi_track;

load('rectangle_dewetting_data64_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');

Nmem2s = hstry.Nmems; Nffi2 = para.Nffi;
pmem2s = hstry.pmem_track;
pffi2s = hstry.pffi_track;

load('rectangle_dewetting_data128_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');

Nmem3s = hstry.Nmems; Nffi3 = para.Nffi;
pmem3s = hstry.pmem_track;
pffi3s = hstry.pffi_track;

load('rectangle_dewetting_data256_cb1e_1_gamma2_0_5_noprod_Ca_0_2.mat');
Nmem4s = hstry.Nmems; Nffi4 = para.Nffi;
pmem4s = hstry.pmem_track;
pffi4s = hstry.pffi_track;

aa = 100; original_dt = 0.02;
order1 = zeros(aa,1); order2 = zeros(aa,1);
for itr = 2:aa+1

ind = itr;
Nmem1 = Nmem1s(ind);
pmem1 = reshape(pmem1s(ind,1:Nmem1,1:2),Nmem1,2);
pffi1 = reshape(pffi1s(ind,1:Nffi1,1:2),Nffi1,2);

ind = (itr-1)*4+1;
Nmem2 = Nmem2s(ind);
pmem2 = reshape(pmem2s(ind,1:Nmem2,1:2),Nmem2,2);
pffi2 = reshape(pffi2s(ind,1:Nffi2,1:2),Nffi2,2);

ind = (itr-1)*16+1;
Nmem3 = Nmem3s(ind);
pmem3 = reshape(pmem3s(ind,1:Nmem3,1:2),Nmem3,2);
pffi3 = reshape(pffi3s(ind,1:Nffi3,1:2),Nffi3,2);

ind = (itr-1)*64+1;
Nmem4 = Nmem4s(ind);
pmem4 = reshape(pmem4s(ind,1:Nmem4,1:2),Nmem4,2);
pffi4 = reshape(pffi4s(ind,1:Nffi4,1:2),Nffi4,2);

error1 = 0;
for i = 1:Nmem1
    dist = abs(p_poly_dist(pmem1(i,:),pmem2));
    error1 = max(error1, dist);
end

for i = 1:Nffi1
    dist = abs(p_poly_dist(pffi1(i,:),pffi2));
    error1 = max(error1, dist);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error2 = 0;
for i = 1:Nmem2
    dist = abs(p_poly_dist(pmem2(i,:),pmem3));
    error2 = max(error2, dist);    
end

for i = 1:Nffi2
    dist = abs(p_poly_dist(pffi2(i,:),pffi3));
    error2 = max(error2, dist);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error3 = 0;
for i = 1:Nmem3
    dist = abs(p_poly_dist(pmem3(i,:),pmem4));
    error3 = max(error3, dist);
end

for i = 1:Nffi3
    dist = abs(p_poly_dist(pffi3(i,:),pffi4));
    error3 = max(error3, dist);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

order1(itr-1) = log(error1/error2)/log(2);

order2(itr-1) = log(error2/error3)/log(2);

if ismember((itr-1)*original_dt,[2])
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