function [memlen] = getlength(p)

N = size(p,1);
diff = p(2:N,1:2)-p(1:N-1,1:2);
memlen = sum(sqrt(diff(:,1).^2+diff(:,2).^2));
