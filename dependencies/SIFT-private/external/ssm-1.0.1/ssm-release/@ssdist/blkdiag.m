function A = blkdiag(varargin)

%@SSDIST/BLKDIAG Block diagonal concatenation of SSDIST objects.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

mat         = cell(1, nargin);
type        = cell(1, nargin);
matf        = cell(1, nargin);
logpf       = cell(1, nargin);
diagmask    = cell(1, nargin);
dmask       = cell(1, nargin);
for i = 1 : nargin
    if isa(varargin{i}, 'ssdist'), d = varargin{i};
    elseif isa(varargin{i}, 'ssmat') || isnumeric(varargin{i}), d = ssdist(varargin{i});
    else error('ssm:ssdist:blkdiag:UnableToConvert', ['Input ' int2str(i) ' cannot be converted to SSDIST class.']);
    end
    mat{i}      = d.ssmat;
    type{i}     = d.type;
    matf{i}     = d.matf;
    logpf{i}    = d.logpf;
    diagmask{i} = d.diagmask;
    dmask{i}    = d.dmask;
end
parent      = blkdiag(mat{:});
A.type      = [type{:}];
A.matf      = [matf{:}];
A.logpf     = [logpf{:}];
A.diagmask  = logical(blkdiag(diagmask{:}));
A.dmask     = [dmask{:}];

%% Register object instance %%
A   = class(A, 'ssdist', parent);


