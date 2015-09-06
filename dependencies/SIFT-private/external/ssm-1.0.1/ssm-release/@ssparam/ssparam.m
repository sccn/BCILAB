function param = ssparam(A, transform, group)

%@SSPARAM/SSPARAM State space model parameter class constructor.
%   param = SSPARAM(w[, transform])
%   param = SSPARAM(name[, transform, group])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin == 0
    %% Default constructor %%
    param.name      = {};
    param.value     = [];
    param.group     = 0;
    param.transform = {};
elseif isa(A, 'ssparam')
    %% Copy constructor %%
    param           = A;
    return;
else
    if isnumeric(A)
        w               = A;
        for i = 1:w, param.name{i} = ['Parameter ' int2str(i)]; end
    elseif ischar(A)
        w               = 1;
        param.name{1}   = A;
    else
        w               = length(A);
        param.name      = A;
    end
    if nargin < 3 || isempty(group), group = [0 w]; else group = [0 cumsum(group)]; end
    if nargin < 2 || isempty(transform), transform = {'identity'}; elseif ischar(transform), temp = cell(1, length(group)-1); [temp{:}] = deal(transform); transform = temp; end
    param.value     = zeros(1, w);
    param.group     = group;
    param.transform = transform;
end

%% Register object instance %%
param   = class(param, 'ssparam');

