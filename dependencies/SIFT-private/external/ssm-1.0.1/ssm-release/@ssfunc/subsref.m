function B = subsref(A, S)

%@SSFUNC/SUBSREF Subscripted reference.
%   .nf returns number of functions.
%   .f returns the functions.
%   .df returns the derivatives of the functions.
%   .horzmask returns the matrix column masks for each function.
%   .vertmask returns the matrix row masks for each function.
%   .fmask returns the variable mask for all functions.
%   .const returns true if SSFUNC (including the underlying SSMAT) is constant.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

switch S(1).type
    case '.'
        switch S(1).subs
            case 'nf',          B = length(A.f);
            case 'f',           B = A.f;
            case 'df',          B = A.df;
            case 'horzmask',    B = A.horzmask;
            case 'vertmask',    B = A.vertmask;
            case 'fmask',       B = A.fmask;
            case 'const',       B = ~any(A.fmask) && A.ssmat.const;
            otherwise,          B = subsref(A.ssmat, S); return;
        end
        if length(S) > 1, B = subsref(B, S(2:end)); end
    otherwise
        B   = subsref(A.ssmat, S);
end



