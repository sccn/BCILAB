function param = remove(param, mask)

%@SSPARAM/REMOVE Remove parameters.
%   param = REMOVE(param, mask)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

remove  = [];
group   = param.group;
for i = 1 : length(param.group)-1
    gmask   = mask(param.group(i)+1:param.group(i+1));
    if any(gmask) && ~all(gmask), mask(param.group(i)+1:param.group(i+1)) = false;
    elseif gmask
        remove          = [remove i];
        group(i+1:end)  = group(i+1:end) - group(i+1) + group(i);
    end
end
group(remove)           = [];
param.name(mask)        = [];
param.value(:, mask)    = [];
param.group             = group;
param.transform(remove) = [];


