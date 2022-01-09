function [D1P1P2] = getStiffonRef1DP1P2

% get (N_i, M_j)_{Ref}
% N_i P1 shape functions for D1 
% M_j P2 shape functions for D2

D1P1P2 = [1/6 1/3 0;
          0 1/3 1/6];