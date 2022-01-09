function [Ls21, Ls22, Ls23, Ls24] = localstiff_movmesh2(p,sp)

x = p(:,1); y = p(:,2);

J = [x(2)-x(1) x(3)-x(1); y(2)-y(1) y(3)-y(1)];
detJ = abs(det(J));
invJ = 1/det(J)*[y(3)-y(1) x(1)-x(3); y(1)-y(2) x(2)-x(1)];
ux = invJ(1,1); uy = invJ(1,2);
vx = invJ(2,1); vy = invJ(2,2);

D2P1der = getStiffonRef2DP1der;

Ls21 = zeros(3); Ls22 = zeros(3); Ls23 = zeros(3); Ls24 = zeros(3);
for i = 1:3
    for j = 1:3
        ref = cell2mat(D2P1der(i,j));
        Ls21(i,j) = sum(sum([ux^2 ux*vx; ux*vx vx^2].*ref));
        Ls22(i,j) = sum(sum([ux*uy ux*vy; uy*vx vx*vy].*ref));
        Ls23(i,j) = sum(sum([uy*ux uy*vx; ux*vy vx*vy].*ref));
        Ls24(i,j) = sum(sum([uy^2 uy*vy; uy*vy vy^2].*ref));
    end
end


Ls21 = sp*detJ*Ls21;
Ls22 = sp*detJ*Ls22;
Ls23 = sp*detJ*Ls23;
Ls24 = sp*detJ*Ls24;
