function model = ssm_binomial(k)

%SSM_BINOMIAL Create SSMODEL object for binomial distribution error model.
%   model = SSM_BINOMIAL(k)
%       k is the number of trials at each time point, or a scalar if the
%           number of trials is stationary.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, error('ssm:ssm_binomial:NotEnoughInputs', 'Insufficient input arguments.'); end
model   = ssmodel(struct('type', 'binomial error', 'k', k), dist_binomial(k), zeros(1, 0), [], [], []);

