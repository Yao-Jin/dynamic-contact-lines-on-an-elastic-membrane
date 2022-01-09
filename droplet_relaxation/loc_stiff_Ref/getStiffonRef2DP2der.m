function [D2P2der] = getStiffonRef2DP2der

% get (Grad N_i, Grad N_j)_{Ref}
% N_i P2 shape functions on the reference triangular

D2P2der = num2cell(zeros(6));

D2P2der(1,1) = {[1/2 1/2; 1/2 1/2]};
D2P2der(1,2) = {[1/6 0; 1/6 0]};
D2P2der(1,3) = {[0 1/6; 0 1/6]};
D2P2der(1,4) = {[-2/3 0; -2/3 0]};
D2P2der(1,5) = {[0 0; 0 0]};
D2P2der(1,6) = {[0 -2/3; 0 -2/3]};

D2P2der(2,2) = {[1/2 0; 0 0]};
D2P2der(2,3) = {[0 -1/6; 0 0]};
D2P2der(2,4) = {[-2/3 -2/3; 0 0]};
D2P2der(2,5) = {[0 2/3; 0 0]};
D2P2der(2,6) = {[0 0; 0 0]};

D2P2der(3,3) = {[0 0; 0 1/2]};
D2P2der(3,4) = {[0 0; 0 0]};
D2P2der(3,5) = {[0 0; 2/3 0]};
D2P2der(3,6) = {[0 0; -2/3 -2/3]};

D2P2der(4,4) = {[4/3 2/3; 2/3 4/3]};
D2P2der(4,5) = {[0 -2/3; -2/3 -4/3]};
D2P2der(4,6) = {[0 2/3; 2/3 0]};

D2P2der(5,5) = {[4/3 2/3; 2/3 4/3]};
D2P2der(5,6) = {[-4/3 -2/3; -2/3 0]};

D2P2der(6,6) = {[4/3 2/3; 2/3 4/3]};

for i = 2:6
    for j = 1:i-1
        D2P2der(i,j) = {(cell2mat(D2P2der(j,i)))'};
    end
end
