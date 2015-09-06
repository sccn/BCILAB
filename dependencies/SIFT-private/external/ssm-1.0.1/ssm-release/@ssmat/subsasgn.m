function A = subsasgn(A, S, B)

%@SSMAT/SUBSASGN Subscripted assignment.
%   .mat sets the stationary part.
%   .dvec sets the dynamic part.
%   (i, j) sets stationary element i, j.
%   (i, j, t) sets element i, j at time t.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

switch S(1).type
    case '.'
        switch S(1).subs
            case 'mat'
                mat     = subsasgn(A.mat, S(2:end), B);
                if ~isequal(size(mat), size(A.mat)), error('ssm:ssmat:subsasgn', 'Invalid subscripted assignment.'); end
                A.mat   = mat;
            case 'dvec'
                dvec    = subsasgn(A.dvec, S(2:end), B);
                if ndims(dvec) ~= ndims(A.dvec) || size(dvec, 1) ~= size(A.dvec, 1), error('ssm:ssmat:subsasgn', 'Invalid subscripted assignment.'); end
                A.dvec  = dvec;
            otherwise, error('ssm:ssmat:subsasgn', 'Invalid subscripted assignment.');
        end
    case '()'
        if length(S(1).subs) >= 3
            if isempty(A.dmmask), error('ssm:ssmat:subsasgn', 'Invalid subscripted assignment.'); end
            T                               = size(A.dvec, 2);
            mat                             = repmat(A.mat, [1 1 T]);
            mat(repmat(A.dmmask, [1 1 T]))  = A.dvec;
            mat2                            = subsasgn(mat, S, B);
            if ~isequal(size(mat2), size(mat)), error('ssm:ssmat:subsasgn', 'Invalid subscripted assignment.'); end
            A.mat                           = mat2(:,:, 1);
            A.mat(A.dmmask)                 = 0;
            A.dvec                          = reshape(mat2(repmat(A.dmmask, [1 1 T])), size(A.dvec));
        else
            mat     = subsasgn(A.mat, S, B);
            if ~isequal(size(mat), size(A.mat)), error('ssm:ssmat:subsasgn', 'Invalid subscripted assignment.'); end
            A.mat   = mat;
        end
    otherwise
        error('ssm:ssmat:subsasgn', 'Invalid subscripted assignment.');
end

