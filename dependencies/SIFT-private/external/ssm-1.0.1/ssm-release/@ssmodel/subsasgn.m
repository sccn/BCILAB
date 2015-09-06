function A = subsasgn(A, S, B)

%@SSMODEL/SUBSASGN Subscripted assignment.
%   .name sets model name.
%   .a1 sets model initial state vector.
%   .P1 sets model initial state variance.
%   .param sets model parameter values.
%   .psi sets transformed model parameter values.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

switch S(1).type
    case '.'
        switch S(1).subs
            case 'name'
                if length(S) == 1
                    if ~ischar(B), error('ssm:ssmodel:subsasgn', 'name must be a string.'); end
                    A.name  = B;
                else A.name = subsasgn(A.name, S(2:end), B);
                end
            case 'c'
                if length(S) == 1
                    if isa(B, 'ssmat') && issta(B) && isequal(size(B), size(A.a1)), A.c = B;
                    elseif isnumeric(B) && isequal(size(B), size(A.a1)), A.c = ssmat(B);
                    else error('ssm:ssmodel:subsasgn', 'c must be a stationary m*1 SSMAT or vector.');
                    end
                else A.c = subsasgn(A.c, S(2:end), B);
                end
            case 'a1'
                if length(S) == 1
                    if isa(B, 'ssmat') && issta(B) && isequal(size(B), size(A.a1)), A.a1 = B;
                    elseif isnumeric(B) && isequal(size(B), size(A.a1)), A.a1 = ssmat(B);
                    else error('ssm:ssmodel:subsasgn', 'a1 must be a stationary m*1 SSMAT or vector.');
                    end
                else A.a1 = subsasgn(A.a1, S(2:end), B);
                end
            case 'P1'
                if length(S) == 1
                    if isa(B, 'ssmat') && issta(B) && isequal(size(B), size(A.P1)), A.P1 = B;
                    elseif isnumeric(B) && isequal(size(B), size(A.P1)), A.P1 = ssmat(B);
                    else error('ssm:ssmodel:subsasgn', 'P1 must be a stationary m*m SSMAT or matrix.');
                    end
                else A.P1 = subsasgn(A.P1, S(2:end), B);
                end
            case 'param'
                if length(S) > 1, B = subsasgn(get(A.psi), S(2:end), B);
                elseif ~isequal(size(B), [1 A.psi.w]), error('ssm:ssmodel:subsasgn', 'Incorrect number of parameters.'); end
                A   = setparam(A, B, false);
            case 'psi'
                if length(S) > 1, B = subsasgn(A.psi.value, S(2:end), B);
                elseif ~isequal(size(B), [1 A.psi.w]), error('ssm:ssmodel:subsasgn', 'Incorrect number of parameters.'); end
                A   = setparam(A, B, true);
            otherwise
                error('ssm:ssmodel:subsasgn', 'Invalid subscripted assignment.');
        end
    otherwise
        error('ssm:ssmodel:subsasgn', 'Invalid subscripted assignment.');
end


