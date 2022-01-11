function [D1P1der] = getStiffonRef1DP1der

% get (Grad M_i, Grad N_j)_{Ref}
% M_i P1 shape function for D1 
% N_j P1 shape function for D1

D1P1der = [1 -1;
           -1 1];