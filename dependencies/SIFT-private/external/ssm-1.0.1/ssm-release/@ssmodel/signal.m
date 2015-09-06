function ycom = signal(alpha, model, t0, varargin)

%@SSMODEL/SIGNAL Retrieve signal components.
%   ycom = SIGNAL(alpha, model[, t0])
%   SIGNAL(..., optname1, optvalue1, optname2, optvalue2, ...)
%       alpha is the state sequence.
%       model is the linear Gaussian model to use.
%       t0 is optional time offset for dynamic Z.
%       ycom is p*n*M where M is number of signal components, unless p == 1,
%           in which case ycom is M*n.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3 || isempty(t0), t0 = 1; end

opt     = setopt([], varargin{:});

checkstate(model, alpha, 'alpha');

if opt.usec
    ycom    = signal_int_c(alpha, model.mcom, getmat_c(model.Z), ~issta(model.Z), t0, true);
else
    n       = size(alpha, 2);
    ncom    = length(model.mcom) - 1;
    Zmat    = getmat(model.Z);
    if issta(model.Z)
        p       = size(Zmat, 1);
        ycom    = zeros(p, n, ncom);
        for i = 1:ncom, ycom(:, :, i) = Zmat(:, model.mcom(i)+1:model.mcom(i+1))*alpha(model.mcom(i)+1:model.mcom(i+1), :); end
    else
        p       = size(Zmat{1}, 1);
        ycom    = zeros(p, n, ncom);
        for t = 1 : n
            Z   = Zmat{t0+t-1};
            for i = 1:ncom, ycom(:, t, i) = Z(:, model.mcom(i)+1:model.mcom(i+1))*alpha(model.mcom(i)+1:model.mcom(i+1), t); end
        end
    end

    if model.p == 1, ycom = permute(ycom, [3 2 1]); end
end
