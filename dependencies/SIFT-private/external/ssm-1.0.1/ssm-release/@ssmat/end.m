function i = end(A, k, n)

%@SSMAT/END Last index.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if n == 1
    i   = numel(A.mat);
else
    if k < 3
        i   = size(A.mat, k);
    elseif k == 3
        i   = size(A.dvec, 2);
    else
        i   = 1;
    end
end

