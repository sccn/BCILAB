function M = getmat(m)

%@SSMAT/GETMAT Get matrix (stationary) or a cell array of matrices (dynamic).

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if isempty(m.dmmask), M = m.mat;
else
    n           = size(m.dvec, 2);
    [M{1:n}]    = deal(m.mat);
    for t = 1 : n, M{t}(m.dmmask) = m.dvec(:, t); end
end

