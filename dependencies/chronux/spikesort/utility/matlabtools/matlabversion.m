function vout = matlabversion
%MATLABVERSION     Sortable version number for current Matlab.
%   V = MATLABVERSION returns a sortable version number for the instance
%   of Matlab invoking MATLABVERSION.  The number is derived by calling  V
%   = ver('Matlab');  V = V.Version;  and converting the result into a
%   number by removing all decimal points after the first.
%
%   MATLABVERSION caches the version number and repeated calls are thus
%   significantly faster than similar calls to VER.
%
%   (Note that the sortable property will break if TMW ever releases an
%   update with more than one digit in the update number.  E.g.,
%   MATLABVERSION assumes that versions such as 7.10.1 do not occur, in
%   keeping with TMW nomenclature to date.)
%
%   See also VER, VERSION.

persistent v

if (isempty(v)),  % only look this up the first time ...
	% Get version string
	v = ver('Matlab');
	v = v.Version;
	
	% Replace all decimal points after the first
	dots = regexp(v,'\.');
	v(dots(2:end)) = '';
	
	% Convert to number
	v = str2num(v);
end

vout = v;