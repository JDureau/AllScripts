function [xm xcov] = UT(Xi,W,C)

[n, kmax] = size(Xi);

xm = zeros(n,1);
for k = 1:kmax
    xm = xm + W(k)*Xi(:,k);
end

xcov = zeros(n,n);

for k = 1:kmax
    xcov = xcov + W(k)*(Xi(:,k)-xm)*((Xi(:,k)-xm)');
end


if nargin == 3
    xcov = xcov + C;
end
