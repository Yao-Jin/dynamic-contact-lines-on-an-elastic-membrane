function [D2P1der] = getStiffonRef2DP1der

% get (Grad N_i, Grad N_j)_{Ref}
% N_i P1 shape functions on the reference triangular

D2P1der = num2cell(zeros(3));

D2P1der(1,1) = {[1/2 1/2; 1/2 1/2]};
D2P1der(1,2) = {[-1/2 0; -1/2 0]};
D2P1der(1,3) = {[0 -1/2; 0 -1/2]};

D2P1der(2,2) = {[1/2 0; 0 0]};
D2P1der(2,3) = {[0 1/2; 0 0]};

D2P1der(3,3) = {[0 0; 0 1/2]};

for i = 2:3
    for j = 1:i-1
        D2P1der(i,j) = {(cell2mat(D2P1der(j,i)))'};
    end
end
