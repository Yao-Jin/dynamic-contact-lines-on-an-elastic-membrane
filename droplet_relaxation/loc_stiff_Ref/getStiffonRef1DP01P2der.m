function [D1P01P2der] = getStiffonRef1DP01P2der

% get (M_i, Grad N_j)_{Ref}
% M_i P0,P1 shape function for D1 
% N_j P2 shape function for D1 

D1P01P2der = [-1 0 1;
              -5/6 2/3 1/6;
              -1/6 -2/3 5/6];