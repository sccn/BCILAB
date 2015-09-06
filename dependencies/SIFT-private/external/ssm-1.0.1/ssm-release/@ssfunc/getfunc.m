function y = getfunc(f, x, t)

%@SSFUNC/GETFUNC Use state space function.
%   y = GETFUNC(f, x, t)
%       x is a m*N matrix.
%       t is a scalar.
%       y is a p*N matrix.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

M   = getmat(f.ssmat);
if ~issta(f.ssmat)
    if getn(f.ssmat) == 1, M = M{1};
    else M = M{t};
    end
end
y   = M*x;

for i = 1 : length(f.f)
    for j = 1 : size(y, 2)
        y(f.vertmask(:, i), j)  = y(f.vertmask(:, i), j) + f.f{i}(x(f.horzmask(:, i), j), t);
    end
end

