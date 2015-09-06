function model = ssm_exp()

%SSM_EXP Create SSMODEL object for exponential distribution error model.
%   model = SSM_EXP()

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

model   = ssmodel('exponential error', dist_exp, zeros(1, 0), [], [], []);

