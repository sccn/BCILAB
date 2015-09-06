function A = horzcat(varargin)

%@SSFUNC/HORZCAT Horizontal concatenation of SSFUNC objects.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

mat         = cell(1, nargin);
f2          = cell(1, nargin);
df          = cell(1, nargin);
horzmask    = cell(1, nargin);
vertmask    = cell(1, nargin);
fmask       = cell(1, nargin);
for i = 1 : nargin
    if isa(varargin{i}, 'ssfunc'), f = varargin{i};
    elseif isa(varargin{i}, 'ssmat') || isnumeric(varargin{i}), f = ssfunc(varargin{i});
    else error('ssm:ssfunc:horzcat:UnableToConvert', ['Input ' int2str(i) ' cannot be converted to SSFUNC class.']);
    end
    mat{i}      = f.ssmat;
    f2{i}       = f.f;
    df{i}       = f.df;
    horzmask{i} = f.horzmask;
    vertmask{i} = f.vertmask;
    fmask{i}    = f.fmask;
end
parent      = [mat{:}];
A.f         = [f2{:}];
A.df        = [df{:}];
A.horzmask  = logical(blkdiag(horzmask{:}));
A.vertmask  = [vertmask{:}];
A.fmask     = [fmask{:}];

%% Register object instance %%
A   = class(A, 'ssfunc', parent);


