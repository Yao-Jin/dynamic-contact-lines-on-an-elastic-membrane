function [D1P01P2] = getStiffonRef1DP01P2

% get (M_i, N_j)_{Ref}
% M_i P0,P1 shape functions for D1
% N_j P2 shape functions for D1 

D1P01P2 = [1/6 2/3 1/6;
           1/6 1/3 0;
           0 1/3 1/6];