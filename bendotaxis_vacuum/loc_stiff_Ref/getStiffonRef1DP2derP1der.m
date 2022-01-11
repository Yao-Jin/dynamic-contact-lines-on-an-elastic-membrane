function [D1P2derP1der] = getStiffonRef1DP2derP1der

% get (partial M_i, partial N_j)_{Ref}
% M_i P2 shape functions for D1
% N_j P1 shape functions for D1

D1P2derP1der = [1 -1; 0 0; -1 1];