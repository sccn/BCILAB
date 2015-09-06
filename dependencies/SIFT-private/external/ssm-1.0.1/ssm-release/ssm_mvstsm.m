function model = ssm_mvstsm(p, cov, lvl, seas, s, cycle, varargin)

%SSM_MVSTSM Create SSMODEL object for multivariate structural time series models.
%   model = SSM_MVSTSM(p, cov, lvl, seas, s[, cycle, x, dep])
%       p is the number of variables.
%       cov specifies complete covariance if true, or complete independence if
%           false, extended to a vector where needed.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 5, error('ssm:ssm_mvstsm:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~ischar(lvl), error('ssm:ssm_mvstsm:InputError', 'lvl must be ''level'' or ''trend''.'); end
if nargin < 6 || isempty(cycle), cycle = false; elseif ~islogical(cycle), error('ssm:ssm_mvstsm:InputError', 'cycle must be logical.'); end

switch lvl
    case 'level', lvl = ssm_mvllm(p, cov(1:2)); k = 3;
    case 'trend', lvl = ssm_mvllt(p, cov(1:3)); k = 4;
    otherwise, lvl = ssm_gaussian(p, cov(1)); k = 2;
end
model   = [lvl ssm_mvseasonal(p, cov(k), seas, s)];
if cycle, model = [model ssm_mvcycle(p, cov(k+1))]; end
if nargin >= 7, model = [model ssm_mvreg(p, varargin{:})]; end

