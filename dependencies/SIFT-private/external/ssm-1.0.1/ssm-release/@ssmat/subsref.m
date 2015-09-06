function B = subsref(A, S)

%@SSMAT/SUBSREF Subscripted reference.
%   .mat returns the stationary part.
%   .mmask returns the variable mask.
%   .n returns time duration.
%   .dmmask returns the dynamic mask.
%   .d returns nnz(dmmask), i.e. number of dynamic elements.
%   .dvec returns the dynamic vector series.
%   .dvmask returns the dynamic vector variable mask.
%   .const returns true if SSMAT is constant.
%   .sta returns true if SSMAT is stationary.
%   (i, j) returns stationary matrix element i, j as a 2-d matrix.
%   (i, j, t) returns matrix element i, j at time t as a 3-d matrix.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

switch S(1).type
    case '.'
        switch S(1).subs
            case 'mat',     B = A.mat;
            case 'mmask',   if isempty(A.mmask), B = false(size(A.mat)); else B = A.mmask; end
            case 'n',       B = getn(A);
            case 'dmmask',  if isempty(A.dmmask), B = false(size(A.mat)); else B = A.dmmask; end
            case 'd',       B = size(A.dvec, 1);
            case 'dvec',    B = A.dvec;
            case 'dvmask',  if isempty(A.dvmask), B = false(size(A.dvec, 1), 1); else B = A.dvmask; end
            case 'const',   B = isconst(A) && isdconst(A);
            case 'sta',     B = issta(A);
            otherwise,      error('ssm:ssmat:subsref', 'Invalid subscripted reference.');
        end
        if length(S) > 1, B = subsref(B, S(2:end)); end
    case '()'
        if length(S(1).subs) >= 3
            if isempty(A.dmmask)
                if isnumeric(S(1).subs{3}), T = max(S(1).subs{3}); else T = 1; end
                B   = repmat(A.mat, [1 1 T]);
            else
                T                               = size(A.dvec, 2);
                B                               = repmat(A.mat, [1 1 T]);
                B(repmat(A.dmmask, [1 1 T]))    = A.dvec;
            end
        else
            B   = A.mat;
        end
        B   = subsref(B, S);
    otherwise
        error('ssm:ssmat:subsref', 'Invalid subscripted reference.');
end

