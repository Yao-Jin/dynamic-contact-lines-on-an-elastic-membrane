function [D2P01P2der] = getStiffonRef2DP01P2der

% get (M_i, Grad N_j)_{Ref}
% M_i P0+P1 shape functions on the reference triangular
% N_j P2 shape functions on the reference triangular

D2P01P2der = num2cell(zeros(4,6));

% (P0, der P2)
D2P01P2der(1,1) = {[-1/6 -1/6]};
D2P01P2der(1,2) = {[1/6 0]};
D2P01P2der(1,3) = {[0 1/6]};
D2P01P2der(1,4) = {[0 -2/3]};
D2P01P2der(1,5) = {[2/3 2/3]};
D2P01P2der(1,6) = {[-2/3 0]};

% (P1, der P2)
D2P01P2der(2,1) = {[-1/6 -1/6]};
D2P01P2der(2,2) = {[0 0]};
D2P01P2der(2,3) = {[0 0]};
D2P01P2der(2,4) = {[1/6 -1/6]};
D2P01P2der(2,5) = {[1/6 1/6]};
D2P01P2der(2,6) = {[-1/6 1/6]};

D2P01P2der(3,1) = {[0 0]};
D2P01P2der(3,2) = {[1/6 0]};
D2P01P2der(3,3) = {[0 0]};
D2P01P2der(3,4) = {[-1/6 -1/3]};
D2P01P2der(3,5) = {[1/6 1/3]};
D2P01P2der(3,6) = {[-1/6 0]};

D2P01P2der(4,1) = {[0 0]};
D2P01P2der(4,2) = {[0 0]};
D2P01P2der(4,3) = {[0 1/6]};
D2P01P2der(4,4) = {[0 -1/6]};
D2P01P2der(4,5) = {[1/3 1/6]};
D2P01P2der(4,6) = {[-1/3 -1/6]};
