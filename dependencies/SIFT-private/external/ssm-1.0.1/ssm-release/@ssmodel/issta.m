function t = issta(model)

%@SSMODEL/ISSTA True for stationary state space models.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

t   = issta(model.H) && issta(model.Z) && issta(model.T) && issta(model.R) && issta(model.Q) && issta(model.c);

