function [Ls11, Ls12, Ls13, Ls14] = localstiff_bulk1(p,eta)

x = p(:,1); y = p(:,2);

J = [x(2)-x(1) x(3)-x(1); y(2)-y(1) y(3)-y(1)];
detJ = abs(det(J));
invJ = 1/det(J)*[y(3)-y(1) x(1)-x(3); y(1)-y(2) x(2)-x(1)];
Xx = invJ(1,1); Xy = invJ(1,2);
Yx = invJ(2,1); Yy = invJ(2,2);


D2P2der = getStiffonRef2DP2der;

Ls11 = zeros(6); Ls12 = zeros(6); Ls13 = zeros(6); Ls14 = zeros(6);
for i = 1:6
    for j = 1:6
        ref = cell2mat(D2P2der(i,j));
        Ls11(i,j) = sum(sum([2*Xx^2+Xy^2 2*Xx*Yx+Xy*Yy; 2*Xx*Yx+Xy*Yy 2*Yx^2+Yy^2].*ref));
        Ls12(i,j) = sum(sum([Xy*Xx Xy*Yx; Xx*Yy Yx*Yy].*ref));
        Ls13(i,j) = sum(sum([Xx*Xy Xx*Yy; Xy*Yx Yx*Yy].*ref));
        Ls14(i,j) = sum(sum([2*Xy^2+Xx^2 2*Xy*Yy+Xx*Yx; 2*Xy*Yy+Xx*Yx 2*Yy^2+Yx^2].*ref));
    end
end

Ls11 = eta*detJ*Ls11;
Ls12 = eta*detJ*Ls12;
Ls13 = eta*detJ*Ls13;
Ls14 = eta*detJ*Ls14;