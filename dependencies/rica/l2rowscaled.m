function [y] = l2rowscaled(x, alpha)

normeps = 1e-5;
epssumsq = sum(x.^2,2) + normeps;   

l2rows=sqrt(epssumsq)*alpha;
y=bsxfun(@rdivide,x,l2rows);
