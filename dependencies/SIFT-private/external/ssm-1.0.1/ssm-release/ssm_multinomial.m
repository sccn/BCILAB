function model = ssm_multinomial(h, k)

%SSM_MULTINOMIAL Create SSMODEL object for multinomial distribution error model.
%   model = SSM_MULTINOMIAL(h, k)
%       h is the number of cells.
%       k is the number of trials at each time point, or a scalar if the
%           number of trials is stationary.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2, error('ssm:ssm_multinomial:NotEnoughInputs', 'Insufficient input arguments.'); end
model   = ssmodel(struct('type', 'multinomial error', 'h', h, 'k', k), dist_multinomial(h, k), zeros(h-1, 0), [], [], []);

