function model = setparam(model, psi, transformed)

%@SSMODEL/SETPARAM Set model parameter values.
%   model = SETPARAM(model, psi[, transformed])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3, transformed = false; end

%% Preprocess parameter values %%
if transformed, model.psi.value = psi;
else
    model.psi       = set(model.psi, psi);
    psi             = model.psi.value;
end

%% Get update adjacency masks %%
Ain             = model.A ~= 0;
Aout            = model.A > 0;
Azero           = ~any(Ain, 1); % We don't need these columns of A
Ain(:, Azero)   = [];
Aout(:, Azero)  = [];
adj_M           = logical(eye(18));
adj_M(:, Azero) = [];
adj_H           = find(adj_M(1, :));
adj_Z           = find(adj_M(2, :));
adj_T           = find(adj_M(3, :));
adj_R           = find(adj_M(4, :));
adj_Q           = find(adj_M(5, :));
adj_c           = find(adj_M(6, :));
adj_a1          = find(adj_M(7, :));
adj_P1          = find(adj_M(8, :));
adj_Hd          = find(adj_M(9, :));
adj_Zd          = find(adj_M(10, :));
adj_Td          = find(adj_M(11, :));
adj_Rd          = find(adj_M(12, :));
adj_Qd          = find(adj_M(13, :));
adj_cd          = find(adj_M(14, :));
adj_Hng         = find(adj_M(15, :));
adj_Qng         = find(adj_M(16, :));
adj_Znl         = find(adj_M(17, :));
adj_Tnl         = find(adj_M(18, :));

%% Get model matrix update values %%
Result  = cell(size(Ain));
for i = 1 : length(model.func), [Result{i, Ain(i, :)}] = model.func{i}(psi(model.pmask{i})); end

%% Update model matrices %%
%%%% TODO: the constant checks could be replaced by checking adj_* instead
if ~isconst(model.H), model.H = setmat(model.H, vertcat(Result{Aout(:, adj_H), adj_H})); end
if ~isconst(model.Z), model.Z = setmat(model.Z, vertcat(Result{Aout(:, adj_Z), adj_Z})); end
if ~isconst(model.T), model.T = setmat(model.T, vertcat(Result{Aout(:, adj_T), adj_T})); end
if ~isconst(model.R), model.R = setmat(model.R, vertcat(Result{Aout(:, adj_R), adj_R})); end
if ~isconst(model.Q), model.Q = setmat(model.Q, vertcat(Result{Aout(:, adj_Q), adj_Q})); end
if ~isconst(model.c), model.c = setmat(model.c, vertcat(Result{Aout(:, adj_c), adj_c})); end
if ~isconst(model.a1), model.a1 = setmat(model.a1, vertcat(Result{Aout(:, adj_a1), adj_a1})); end
if ~isconst(model.P1), model.P1 = setmat(model.P1, vertcat(Result{Aout(:, adj_P1), adj_P1})); end
if ~isdconst(model.H), model.H = setdvec(model.H, vertcat(Result{Aout(:, adj_Hd), adj_Hd})); end
if ~isdconst(model.Z), model.Z = setdvec(model.Z, vertcat(Result{Aout(:, adj_Zd), adj_Zd})); end
if ~isdconst(model.T), model.T = setdvec(model.T, vertcat(Result{Aout(:, adj_Td), adj_Td})); end
if ~isdconst(model.R), model.R = setdvec(model.R, vertcat(Result{Aout(:, adj_Rd), adj_Rd})); end
if ~isdconst(model.Q), model.Q = setdvec(model.Q, vertcat(Result{Aout(:, adj_Qd), adj_Qd})); end
if ~isdconst(model.c), model.c = setdvec(model.c, vertcat(Result{Aout(:, adj_cd), adj_cd})); end
if isa(model.H, 'ssdist') && ~isdistconst(model.H), model.H = setdist(model.H, vertcat(Result{Aout(:, adj_Hng), adj_Hng})); end
if isa(model.Q, 'ssdist') && ~isdistconst(model.Q), model.Q = setdist(model.Q, vertcat(Result{Aout(:, adj_Qng), adj_Qng})); end
if isa(model.Z, 'ssfunc') && ~isfconst(model.Z), model.Z = setfunc(model.Z, vertcat(Result{Aout(:, adj_Znl), adj_Znl})); end
if isa(model.T, 'ssfunc') && ~isfconst(model.T), model.T = setfunc(model.T, vertcat(Result{Aout(:, adj_Tnl), adj_Tnl})); end

