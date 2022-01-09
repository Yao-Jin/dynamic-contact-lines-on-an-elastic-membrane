function [D1P2P2] = getStiffonRef1DP2P2

% get (N_i, N_j)_{Ref}
% N_{i,j} P2 shape functions for D1 

D1P2P2 = [2/15 1/15 -1/30;
        1/15 8/15 1/15;
       -1/30 1/15 2/15];