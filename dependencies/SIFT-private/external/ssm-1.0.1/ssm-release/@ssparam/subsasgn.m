function A = subsasgn(A, S, B)

%@SSPARAM/SUBSASGN Subscripted assignment.
%   .value sets the transformed parameter values.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

switch S(1).type
    case '.'
        switch S(1).subs
            case 'value'
                if length(S) == 1, value = B;
                else value = subsasgn(A.value, S(2:end), B);
                end
                if ~isequal(size(value), size(A.value)), error('ssm:ssparam:subsasgn', 'Invalid subscripted assignment.'); end
                A.value     = value;
            otherwise
                error('ssm:ssparam:subsasgn', 'Invalid subscripted assignment.');
        end
    otherwise
        error('ssm:ssparam:subsasgn', 'Invalid subscripted assignment.');
end




