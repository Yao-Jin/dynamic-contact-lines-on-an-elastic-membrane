function [D2P2P2] = getStiffonRef2DP2P2

% get (M_i, Grad N_j)_{Ref}
% M_i P0+P1 shape functions on the reference triangular
% N_j P2 shape functions on the reference triangular

D2P2P2 = [1/60 -1/360 -1/360 0 -1/90 0;
         -1/360 1/60 -1/360 0 0 -1/90;
         -1/360 -1/360 1/60 -1/90 0 0;
         0 0 -1/90 4/45 2/45 2/45;
         -1/90 0 0 2/45 4/45 2/45;
         0 -1/90 0 2/45 2/45 4/45];
