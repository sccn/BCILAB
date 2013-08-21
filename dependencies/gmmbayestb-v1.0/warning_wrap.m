% WARNING_WRAP()  warning function wrapper
% to allow Matlab R13 style warning calls in Matlab R12

% $Name:  $
% $Id: warning_wrap.m,v 1.1 2004/11/02 08:32:22 paalanen Exp $

% Pekka Paalanen, 2004

function [] = warning_wrap(varargin);

old_version = strcmp(version('-release'), '12');

if old_version
	if nargin > 1
		warning(varargin{2});
	else
		warning(varargin{1});
	end
else
	% Assume version Matlab R13
	warning(varargin{:});
end
