function [D1P1P2der] = getStiffonRef1DP1P2der

% get (M_i, Grad N_j)_{Ref}
% M_i P1 shape function for D1 
% N_j P2 shape function for D1 

D1P1P2der = [-5/6 2/3 1/6;
             -1/6 -2/3 5/6];