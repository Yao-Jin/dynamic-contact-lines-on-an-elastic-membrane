function [D2P2P1] = getStiffonRef2DP2P1

% get (M_i, Grad N_j)_{Ref}
% M_i P0+P1 shape functions on the reference triangular
% N_j P2 shape functions on the reference triangular

D2P2P1 = [1/60 -1/120 -1/120;
          -1/120  1/60 -1/120;
          -1/120 -1/120 1/60;
          1/15 1/15 1/30;
          1/30 1/15 1/15;
          1/15 1/30 1/15];
