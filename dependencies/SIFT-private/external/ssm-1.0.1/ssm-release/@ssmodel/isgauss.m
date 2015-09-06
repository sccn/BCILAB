function t = isgauss(model)

%@SSMODEL/ISGAUSS True for Gaussian state space models.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

t   = ~isa(model.H, 'ssdist') && ~isa(model.Q, 'ssdist');

