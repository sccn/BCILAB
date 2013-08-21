function x=MakeRDSquares(N,nbs,Dyna)
% NESTA Version 1.1


if nargin <3
    Dyna = 40;
end

if nargin <2
    nbs = 5;
end

if nargin <1
    N=256;
end

lmin = floor(8);
lmax = floor(N/4);
x = zeros(N,N);

for ll=1:nbs
    ndx = 1+floor((N-lmax-1)*rand(1,1));
    lx = min(N-ndx-1,floor(lmin + (lmax-lmin)*rand(1,1)));
    ndy = 1+floor((N-lmax-1)*rand(1,1));
    ly = min(N-ndy-1,floor(lmin + (lmax-lmin)*rand(1,1)));
    x(ndx:ndx+lx-1,ndy:ndy+ly-1) = 1+10^(Dyna/20)*rand(1,1);
end
ind = find(x > 0.5);
x(ind) = x(ind)-min(x(ind));
x(ind) = x(ind)/max(x(ind))*(10^(Dyna/20)-1)+1;
 
