function [D1P2P1] = getStiffonRef1DP2P1

% get (N_i, M_j)_{Ref}
% N_i P2 shape functions for D1 
% M_j P1 shape functions for D1 

D1P2P1 = [1/6 0; 1/3 1/3; 0 1/6];