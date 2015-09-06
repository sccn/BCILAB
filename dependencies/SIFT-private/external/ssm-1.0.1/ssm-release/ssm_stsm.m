function model = ssm_stsm(lvl, seas, s, cycle, varargin)

%SSM_STSM Create SSMODEL object for structural time series models.
%   model = SSM_STSM(lvl, seas, s[, cycle, x, varname])
%       lvl is 'level' or 'trend'.
%       seas is the seasonal type (see ssm_seasonal).
%       s is the seasonal period.
%       Set cycle to true if there is a cycle component in the model.
%       x is explanatory variables.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3, error('ssm:ssm_stsm:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~ischar(lvl), error('ssm:ssm_stsm:InputError', 'lvl must be ''level'' or ''trend''.'); end
if nargin < 4, cycle = false; end

switch lvl
    case 'level', model = ssm_llm;
    case 'trend', model = ssm_llt;
    otherwise, model = ssm_gaussian;
end
model   = [model ssm_seasonal(seas, s)];
if cycle, model = [model ssm_cycle]; end
if nargin >= 5, model = [model ssm_reg(varargin{:})]; end

