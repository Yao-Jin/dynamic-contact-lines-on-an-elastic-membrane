function [Ls21, Ls22] = localstiff_bulk2(p)

x = p(:,1); y = p(:,2);

J = [x(2)-x(1) x(3)-x(1); y(2)-y(1) y(3)-y(1)];
detJ = abs(det(J));
invJ = 1/det(J)*[y(3)-y(1) x(1)-x(3); y(1)-y(2) x(2)-x(1)];
ux = invJ(1,1); uy = invJ(1,2);
vx = invJ(2,1); vy = invJ(2,2);

D2P01P2der = getStiffonRef2DP01P2der;

Ls21 = zeros(4,6); Ls22 = zeros(4,6);
for i = 1:4
    for j = 1:6
        ref = cell2mat(D2P01P2der(i,j));
        Ls21(i,j) = sum([ux vx].*ref);
        Ls22(i,j) = sum([uy vy].*ref);
    end
end

Ls21 = detJ*Ls21;
Ls22 = detJ*Ls22;


