function [param pmask] = horzcat(varargin)

%@SSPARAM/HORZCAT Horizontal concatenation of SSPARAM objects.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

N           = 1;
name        = cell(1, nargin);
value       = cell(1, nargin);
group       = 0;
transform   = cell(1, nargin);
fmask       = false(1, 0);
pmask       = cell(1, nargin);
for i = 1 : nargin
    p               = varargin{i}; if ~isa(p, 'ssparam'), error('ssm:ssparam_horzcat:UnableToConvert', ['Input ' int2str(i) ' cannot be converted to SSPARAM class.']); end
    N               = max([N size(p.value, 1)]);
    name{i}         = p.name;
    value{i}        = p.value;
    group           = [group group(end) + p.group(2:end)];
    transform{i}    = p.transform;
    pmask{i}        = [fmask true(1, length(p.name))];
    fmask           = [fmask false(1, length(p.name))];
end
if N > 1, for i = 1:nargin, value{i} = [value{i}; zeros(N-size(value{i}, 1), size(value{i}, 2))]; end, end
param.name      = [name{:}];
param.value     = [value{:}];
param.group     = group;
param.transform = [transform{:}];
if nargout > 1
    w   = length(param.name);
    for i = 1 : nargin
        pmask{i}    = [pmask{i} false(1, w-length(pmask{i}))];
    end
end

%% Register object instance %%
param   = class(param, 'ssparam');




