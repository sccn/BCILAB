function t = islinear(model)

%@SSMODEL/ISLINEAR True for linear state space models.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

t   = ~isa(model.Z, 'ssfunc') && ~isa(model.T, 'ssfunc');

