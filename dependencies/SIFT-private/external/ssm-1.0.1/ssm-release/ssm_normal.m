function model = ssm_normal(varargin)

%SSM_NORMAL Create SSMODEL object for Gaussian noise model.
%   model = SSM_NORMAL(p[, cov])
%       p is the number of variables.
%       cov specifies complete covariance if true, or complete independence if
%           false.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

model   = ssm_gaussian(varargin{:});

