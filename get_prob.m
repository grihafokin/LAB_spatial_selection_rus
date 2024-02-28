function [prb, dA, prbR] = get_prob(x, y, stdC, prbA)

xA = (-10*stdC:0.1:10*stdC) + x;
yA = (-10*stdC:0.1:10*stdC) + y;

[X,Y] = meshgrid(xA, yA);
p = mvnpdf([X(:) Y(:)], [x, y], diag([stdC, stdC].^2));
p = reshape(p, size(X));

dA = 0:0.5:4*stdC;
N = length(dA);
prb = zeros(1, N);
for i=1:N
    prb(i) = sum(p((sqrt((X - x).^2+(Y - y).^2) <= dA(i))))*0.1^2;
end

M = length(prbA);
prbR = zeros(1, M);
for i=1:M
    [~, ind] = min(abs(prb - prbA(i)));
    prbR(i) = dA(ind);
end

end