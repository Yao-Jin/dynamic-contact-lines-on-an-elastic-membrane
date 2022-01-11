function [D1P2P2der] = getStiffonRef1DP2P2der

% get (N_i, Grad N_j)_{Ref}
% N_{i,j} P2 shape functions for D1 

D1P2P2der = [-1/2 2/3 -1/6;
           -2/3  0   2/3;
            1/6 -2/3 1/2];