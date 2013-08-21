function [picks,P] = selectPF(tol,n)

halfn = round(n/2);
[I,J] = meshgrid(1:n,1:n);
P = abs(I - halfn*(rand(n) + .5)) < tol & ...
    abs(J - halfn*(rand(n) + .5)) < tol;

P(halfn+1:n,:) = 0;
P(halfn:halfn+1,:) = 1;
P(:,halfn:halfn+1) = 1;
P = ifftshift(P);
P(1,1) = 1;

picks = find(P);