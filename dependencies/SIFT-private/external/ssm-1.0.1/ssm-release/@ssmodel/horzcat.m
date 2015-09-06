function model = horzcat(varargin)

%@SSMODEL/HORZCAT Horizontal concatenation of SSMODEL objects.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

%% Adjacency logical matrix entries %%
adj_M   = logical(eye(18));
adj_H   = adj_M(1, :);
adj_Hd  = adj_M(9, :);
adj_Hng = adj_M(15, :);
adj_allH    = adj_H | adj_Hng | adj_Hd;

%% Combine models %%
p       = -1;
name    = '';
info    = cell(1, nargin);
mcom    = 0;
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
    m           = varargin{i}; if ~isa(m, 'ssmodel'), error('ssm:ssmodel:horzcat:UnableToConvert', ['Input ' int2str(i) ' cannot be converted to SSMODEL class.']); end
    if p < 0, p = m.p; elseif p ~= m.p, error('ssm:ssmodel:horzcat:IncompatibleModel', 'All models must have the same observation vector length.'); end
    if ~isempty(name), if ~isempty(m.name), name = [name ' + ' m.name]; end, else name = m.name; end
    info{i}     = m.info;
    if ~isscalar(m.mcom), mcom = [mcom mcom(end) + m.mcom(2:end)]; end
    Z{i}        = m.Z;
    T{i}        = m.T;
    R{i}        = m.R;
    Q{i}        = m.Q;
    c{i}        = m.c;
    a1{i}       = m.a1;
    P1{i}       = m.P1;
    if i > 1
        m.A(repmat(adj_allH, size(m.A, 1), 1) & m.A == 1)   = -1;
        removal                                             = all(m.A <= 0, 2);
        if any(removal)
            m.A(removal, :)     = [];
            m.func(removal)     = [];
            m.grad(removal)     = [];
            premoval            = sum(vertcat(m.pmask{removal}), 1) > 0; % parameters used by removed functions
            if any(~removal), premoval = premoval & ~(sum(vertcat(m.pmask{~removal}), 1) > 0); end % parameters used by removed functions and not used by any other functions
            m.psi               = remove(m.psi, premoval); % remove said parameters
            m.pmask(removal)    = []; % remove the parameter masks used by removed functions
            for j = 1 : length(m.pmask), m.pmask{j}(premoval) = []; end
        end
    end
    A{i}        = m.A;
    func{i}     = m.func;
    grad{i}     = m.grad;
    psi{i}      = m.psi;
    pmask{i}    = m.pmask;
end
model.name          = name;
model.Hinfo         = varargin{1}.Hinfo;
model.info          = [info{:}];
model.mcom          = mcom;
model.H             = varargin{1}.H;
model.Z             = [Z{:}];
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


