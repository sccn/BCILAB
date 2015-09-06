function valid = shouldbevalid(A)

%@SSPARAM/SHOULDBEVALID State space model parameter should be valid.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

valid_type  = iscell(A.name) && ndims(A.name) == 2 && (isempty(A.name) || size(A.name, 1) == 1) ...
    && isnumeric(A.value) && size(A.name, 2) == size(A.value, 2) ...
    && isnumeric(A.group) && size(A.group, 1) == 1 && A.group(1) == 0 && A.group(end) == length(A.name) && all(diff(A.group) > 0) ...
    && iscell(A.transform) && (isempty(A.transform) || isequal(size(A.transform), size(A.group)-[0 1]));

valid_content   = true;
for i = 1 : length(A.name)
    valid_content = valid_content && ischar(A.name{i});
end
for i = 1 : length(A.transform)
    valid_content = valid_content && ischar(A.transform{i});
end

valid   = valid_type && valid_content;

if ~valid, warning('ssm:ssparam:shouldbevalid', 'Invalid state space model parameter.'); end

