function [D1P2der] = getStiffonRef1DP2der

% get (Grad N_i, Grad N_j)_{Ref}
% N_{i,j} P2 shape functions for D1 

D1P2der = [7/3 -8/3 1/3;
           -8/3 16/3 -8/3;
           1/3 -8/3 7/3];