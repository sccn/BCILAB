function model = ssm_expfamily(b, d2b, id2bdb, c)

%SSM_EXPFAMILY Create SSMODEL object for general exponential family distribution error model.
%   model = SSM_EXPFAMILY(b, d2b, id2bdb, c)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 4, error('ssm:ssm_expfamily:NotEnoughInputs', 'Insufficient input arguments.'); end
model   = ssmodel('exponential family error', dist_expfamily(b, d2b, id2bdb, c), zeros(1, 0), [], [], []);

