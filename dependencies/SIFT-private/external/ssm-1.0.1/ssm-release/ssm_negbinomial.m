function model = ssm_negbinomial(k)

%SSM_NEGBINOMIAL Create SSMODEL object for negative binomial distribution error model.
%   model = SSM_NEGBINOMIAL(k)
%       k is the number of trials at each time point, or a scalar if the
%           number of trials is stationary.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, error('ssm:ssm_negbinomial:NotEnoughInputs', 'Insufficient input arguments.'); end
model   = ssmodel(struct('type', 'negative binomial error', 'k', k), dist_negbinomial(k), zeros(1, 0), [], [], []);

