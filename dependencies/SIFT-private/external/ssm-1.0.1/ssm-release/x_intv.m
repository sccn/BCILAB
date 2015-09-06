function x = x_intv(n, type, tau)

%X_INTV Create regression variables for intervention components.
%   x = X_INTV(n, type, tau)
%       n is the time series length.
%       type specifies intervention type: 'step', 'pulse', 'slope' or 'null'.
%       tau is the intervention onset time.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

switch type
    case 'step'
        x   = [zeros(1, tau-1) ones(1, n-tau+1)];
    case 'pulse'
        x   = [zeros(1, tau-1) 1 zeros(1, n-tau)];
    case 'slope'
        x   = [zeros(1, tau-1) 1:n-tau+1];
    case 'null'
        x   = zeros(1, n);
end

