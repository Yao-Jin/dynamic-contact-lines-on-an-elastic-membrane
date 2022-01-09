function [D1P2P1der] = getStiffonRef1DP2P1der

% get (M_i, partial N_j)_{Ref}
% M_i P2 shape functions for D1
% N_j P1 shape functions for D1

D1P2P1der = [-1/6 1/6; -2/3 2/3; -1/6 1/6];