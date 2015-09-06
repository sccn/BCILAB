function varargout = size(A, varargin)

%@SSMAT/SIZE Size of stationary part of state space matrix.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

[varargout{1:nargout}] = size(A.mat, varargin{:});

