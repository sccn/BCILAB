function B = subsref(A, S)

%@SSPARAM/SUBSREF Subscripted reference.
%   .w returns number of parameters.
%   .N returns number of sets of parameters.
%   .con returns true if parameter values are constrained.
%   .name returns cell array of parameter names.
%   .value returns transformed parameter values.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

switch S(1).type
    case '.'
        switch S(1).subs
            case 'w',       B = A.group(end);
            case 'N',       B = size(A.value, 1);
            case 'con'
                B   = false;
                for i = 1 : length(A.group)-1
                    if any(strcmp(A.transform{i}, {'covariance' 'ar>=3' 'ma>=3'}))
                        B = true;
                        break;
                    end
                end
            case 'name',    B = A.name;
            case 'value',   B = A.value;
            otherwise,      error('ssm:ssparam:subsref', 'Invalid subscripted reference.');
        end
    otherwise
        error('ssm:ssparam:subsref', 'Invalid subscripted reference.');
end
if length(S) > 1, B = subsref(B, S(2:end)); end

