% close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error

site = 'rectangle_wetting_data128_cb1e_1_gamma1_0_5_reprod_Ca_0_2.mat';

load(site);

sln.nu = para.Ca*sln.bargamma;
sln.nu_old = para.Ca*sln.bargamma_old;
hstry.nu_track = para.Ca*hstry.bargamma_track;
% sln = rmfield(sln, 'bargamma');
% sln = rmfield(sln, 'bargamma_old');
% hstry = rmfield(hstry, 'bargamma_track');

save(site);