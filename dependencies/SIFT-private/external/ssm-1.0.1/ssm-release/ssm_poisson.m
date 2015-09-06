function model = ssm_poisson()

%SSM_POISSON Create SSMODEL object for Poisson distribution error model.
%   model = SSM_POISSON()

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

model   = ssmodel('poisson error', dist_poisson, zeros(1, 0), [], [], []);

