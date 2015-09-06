function B = subsref(A, S)

%@SSDIST/SUBSREF Subscripted reference.
%   .nd returns the number of distributions.
%   .type returns the type of each distribution.
%   .matf returns functions to calculate Gaussian approximating matrices for each distribution.
%   .logpf returns functions to calculate log probabilities for each distribution.
%   .diagmask returns the matrix diagonal masks for each distribution.
%   .dmask returns variable mask for all distributions.
%   .const returns true if SSDIST (including the underlying SSMAT) is constant.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

switch S(1).type
    case '.'
        switch S(1).subs
            case 'nd',          B = length(A.type);
            case 'type',        B = A.type;
            case 'matf',        B = A.matf;
            case 'logpf',       B = A.logpf;
            case 'diagmask',    B = A.diagmask;
            case 'dmask',       B = A.dmask;
            case 'const',       B = ~any(A.dmask) && A.ssmat.const;
            otherwise,          B = subsref(A.ssmat, S); return;
        end
        if length(S) > 1, B = subsref(B, S(2:end)); end
    otherwise
        B   = subsref(A.ssmat, S);
end



