function [Z Zdmmask Zdvec T R] = mat_mvintv(p, n, type, tau)

%MAT_MVINTV Create base matrices for multivariate intervention component.
%   [Z Zdmmask Zdvec T R] = MAT_MVINTV(p, n, type, tau)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

mask                            = logical(eye(p));
mask(:, strcmp(type, 'null'))   = [];
Zdvec   = cell(1, p);
for i = 1 : p
    switch type{i}
        case 'step', Zdvec{i} = [zeros(1, tau(i)-1) ones(1, n-tau(i)+1)];
        case 'pulse', Zdvec{i} = [zeros(1, tau(i)-1) 1 zeros(1, n-tau(i))];
        case 'slope', Zdvec{i} = [zeros(1, tau(i)-1) 1:n-tau(i)+1];
        case 'null', Zdvec{i} = zeros(0, n);
    end
end
m       = size(mask, 2);
Z       = zeros(size(mask));
Zdmmask = mask;
Zdvec   = vertcat(Zdvec{:});
T       = eye(m);
R       = zeros(m, 0);

