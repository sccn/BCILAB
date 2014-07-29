function [z,zxind,zyind,xflag,yflag] = cadaunion(x,y,m,n)
% CADA derivative UNION
% 
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

xrows = x(:,1); xcols = x(:,2); nzx = size(x,1);
yrows = y(:,1); ycols = y(:,2); nzy = size(y,1);

x = sparse(xrows,xcols,-ones(nzx,1),m,n);
y = sparse(yrows,ycols,2*ones(nzy,1),m,n);

z = x+y;
% Get nonzero locs of z deriv
[zrows,zcols,zind] = find(z);
if size(zrows,2) > 1
    zrows = zrows'; zcols = zcols'; zind = zind';
end
z = [zrows,zcols]; nzz = length(zrows);

% Get which indices of z correspond to x
zxind = find(zind<2);
if nzz == nzx; xflag = 1; else xflag = 0; end

% Get which indices of z correspond to y
zyind = find(zind>0);
if nzz == nzy; yflag = 1; else yflag = 0; end