function d = setdist(d, distf)

%@SSDIST/SETDIST Update state space distribution.
%   d = SETDIST(d, distf)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

[d.matf{d.dmask}]  = distf{:, 1};
[d.logpf{d.dmask}] = distf{:, 2};

