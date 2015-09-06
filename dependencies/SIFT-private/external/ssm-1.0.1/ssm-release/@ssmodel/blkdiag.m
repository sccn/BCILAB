function model = blkdiag(varargin)

%@SSMODEL/BLKDIAG Block diagonal concatenation of SSMODEL objects.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

%% Combine independent models %%
name    = '';
Hinfo   = struct('type', 'Multiple noise', 'noise', {cell(1, nargin)});
info    = cell(1, nargin);
mcom    = 0;
H       = cell(1, nargin);
Z       = cell(1, nargin);
T       = cell(1, nargin);
R       = cell(1, nargin);
Q       = cell(1, nargin);
c       = cell(1, nargin);
a1      = cell(1, nargin);
P1      = cell(1, nargin);
A       = cell(1, nargin);
func    = cell(1, nargin);
grad    = cell(1, nargin);
psi     = cell(1, nargin);
pmask   = cell(1, nargin);
for i = 1 : nargin
    m           = varargin{i}; if ~isa(m, 'ssmodel'), error('ssm:ssmodel:blkdiag:UnableToConvert', ['Input ' int2str(i) ' cannot be converted to SSMODEL class.']); end
    if ~isempty(name), if ~isempty(m.name), name = [name ' * ' m.name]; end, else name = m.name; end
    Hinfo.noise{i}  = m.Hinfo;
    info{i}         = m.info;
    if ~isscalar(m.mcom), mcom = [mcom mcom(end) + m.mcom(2:end)]; end
    H{i}        = m.H;
    Z{i}        = m.Z;
    T{i}        = m.T;
    R{i}        = m.R;
    Q{i}        = m.Q;
    c{i}        = m.c;
    a1{i}       = m.a1;
    P1{i}       = m.P1;
    A{i}        = m.A;
    func{i}     = m.func;
    grad{i}     = m.grad;
    psi{i}      = m.psi;
    pmask{i}    = m.pmask;
end
model.name          = name;
model.Hinfo         = Hinfo;
model.info          = [info{:}];
model.mcom          = mcom;
model.H             = blkdiag(H{:});
model.Z             = blkdiag(Z{:});
model.T             = blkdiag(T{:});
model.R             = blkdiag(R{:});
model.Q             = blkdiag(Q{:});
model.c             = vertcat(c{:});
model.a1            = vertcat(a1{:});
model.P1            = blkdiag(P1{:});
model.A             = vertcat(A{:});
model.func          = [func{:}];
model.grad          = [grad{:}];
[model.psi pmask2]  = horzcat(psi{:});
for i = 1 : nargin
    for j = 1 : length(pmask{i})
        temp        = pmask2{i};
        temp(temp)  = pmask{i}{j};
        pmask{i}{j} = temp;
    end
end
model.pmask         = [pmask{:}];
model.p             = size(model.H, 1);
model.n             = max([model.H.n model.Z.n model.T.n model.R.n model.Q.n model.c.n]);

%% Register object instance %%
model   = class(model, 'ssmodel');


