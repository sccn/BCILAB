function model = ssm_arimacom(d, D, s, phi, theta, ksivar)

%SSM_ARIMACOM Create SSMODEL object for ARIMA component model.
%   model = SSM_ARIMACOM(d, D, s, phi, theta, ksivar)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 6, error('ssm:ssm_arimacom:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isnumeric(d) || ~isscalar(d), error('ssm:ssm_arimacom:InputError', 'd must be scalar.'); end
if ~isnumeric(D) || ~isscalar(D), error('ssm:ssm_arimacom:InputError', 'D must be scalar.'); end
if ~isscalar(s) || ~isnumeric(s) || s < 1, error('ssm:ssm_arimacom:InputError', 's must be a scalar no less than one.'); end
if isempty(phi), phi = {}; elseif isnumeric(phi), phi = {phi}; elseif ~iscell(phi), error('ssm:ssm_arimacom:InputError', 'phi must be a cell array of vectors.'); end
if ~iscell(theta), error('ssm:ssm_arimacom:InputError', 'theta must be a cell array of vectors.'); end
if ~isnumeric(ksivar), error('ssm:ssm_arimacom:InputError', 'ksivar must be numeric.'); end
M       = length(phi);
param   = ksivar(end);
dD      = 0;
model   = ssm_gaussian;
% Ordinary difference part
if d > 0
    dD      = dD + 1;
    q       = length(theta{dD});
    model   = [model ssm_arima(0, d, q, false)];
    param   = [param theta{dD} ksivar(dD)];
end
% Seasonal difference part
if D > 0
    dD      = dD + 1;
    q       = length(theta{dD});
    model   = [model ssm_sumarma(0, q, D, s, false)];
    param   = [param theta{dD} ksivar(dD)];
end
% Autoregressive parts
for i = 1 : M
    p       = length(phi{i});
    q       = length(theta{i+dD});
    model   = [model ssm_arma(p, q, false)];
    param   = [param phi{i} theta{i+dD} ksivar(i+dD)];
end
% Moving average part
if length(theta) > M+dD
    q       = length(theta{end});
    model   = [model ssm_arma(0, q, false)];
    param   = [param theta{end} ksivar(end-1)];
end
% Construct model and set parameters
model       = setparam(model, param, false);

