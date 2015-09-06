%CPHELP		show CPRINTF documentation in a ML web browser
%
%		see also: cprintf
%
%SYNTAX
%		CPHELP;
%
%NOTE
%		the published M-file cphelp.html must be
%		- in the current folder
%		- in a folder in ML's search path
%		use a fixed font for HTML proportional text (preferences)
%		  for proper formatting of the output

% created:
%	us	20-Apr-2009 us@neurol.unizh.ch
% modified:
%	us	11-Jun-2009 08:58:55
%
% localid:	us@USZ|ws-nos-36362|x86|Windows XP|7.8.0.347.R2009a

%-------------------------------------------------------------------------------
function	cphelp(varargin)

		url='cphelp.html';

	if	exist(url,'file')
		web(url);
	else
		disp(sprintf('CPRINTF> help file <%s> not found',url));
	end
end
%-------------------------------------------------------------------------------