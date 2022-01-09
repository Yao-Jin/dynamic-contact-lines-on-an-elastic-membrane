function [Ls11, Ls12, Ls13, Ls14] = localstiff_movmesh1(p,sp)

x = p(:,1); y = p(:,2);

J = [x(2)-x(1) x(3)-x(1); y(2)-y(1) y(3)-y(1)];
detJ = abs(det(J));
invJ = 1/det(J)*[y(3)-y(1) x(1)-x(3); y(1)-y(2) x(2)-x(1)];
ux = invJ(1,1); uy = invJ(1,2);
vx = invJ(2,1); vy = invJ(2,2);


D2P1der = getStiffonRef2DP1der;

Ls11 = zeros(3); Ls12 = zeros(3); Ls13 = zeros(3); Ls14 = zeros(3);
for i = 1:3
    for j = 1:3
        ref = cell2mat(D2P1der(i,j));
        Ls11(i,j) = sum(sum([2*ux^2+uy^2 2*ux*vx+uy*vy; 2*ux*vx+uy*vy 2*vx^2+vy^2].*ref));
        Ls12(i,j) = sum(sum([uy*ux uy*vx; ux*vy vx*vy].*ref));
        Ls13(i,j) = sum(sum([ux*uy ux*vy; uy*vx vx*vy].*ref));
        Ls14(i,j) = sum(sum([2*uy^2+ux^2 2*uy*vy+ux*vx; 2*uy*vy+ux*vx 2*vy^2+vx^2].*ref));
    end
end

Ls11 = sp*detJ*Ls11;
Ls12 = sp*detJ*Ls12;
Ls13 = sp*detJ*Ls13;
Ls14 = sp*detJ*Ls14;
