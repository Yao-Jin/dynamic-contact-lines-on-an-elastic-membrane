close all;
clear;
clc;

% gamma1 = 1; gamma2 = 0.5; gamma3 = 1;
exacttheta12 = pi-acos(0.25)
exacttheta23 = exacttheta12
exacttheta31 = pi-2*(asin(0.25))

% % gamma1 =0.5; gamma2 = 1; gamma3 = 1;
exacttheta12 = pi-acos(0.25)
exacttheta23 = pi-2*(asin(0.25))
exacttheta31 = exacttheta12

% dewetting 1// wetting 2
ind = 2;

kmem = [2.05 1.638];
R1 = 1/kmem(ind);

kffi = [2.025 0.852];
R3 = 1/kffi(ind);

d = [0.4796 0.6012];
num_theta1 = asin(d(ind)/R1);
num_theta3 = asin(d(ind)/R3);
% theta3 = asin(d/R3);
num_theta12 = pi-num_theta1
num_theta23 = pi-num_theta3
num_theta31 = (num_theta1+num_theta3)


