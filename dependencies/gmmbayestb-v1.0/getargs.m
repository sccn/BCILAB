%GETARGS  parse variable argument list into a struct
%
% S = GETARGS(defaultS, varglist)
%
% varglist - a cell array of name, value pairs
% defaultS - struct containing the default values
%
% Example:
%   function foo(par1, varargin);
%   args = struct( 'param1', 0, 'param2', eye(3) );
%   args = getargs( args, varargin );
%   disp(args.param1);
%
%  foo(2, 'param1', 14) will print 14
%
% Author:
%   Pekka Paalanen <paalanen@lut.fi>
%
% $name$
% $Id: getargs.m,v 1.1 2004/11/02 08:32:22 paalanen Exp $

function S = getargs(defaultS, varglist);

if mod(length(varglist),2) ~=0
	error('Odd number of variable parameters');
end

S = defaultS;
i=1;
while i <= length(varglist)
	if isfield(S, varglist{i})
		% for Matlab R12
		%S = setfield(S, varglist{i}, varglist{i+1});
		
		% for Matlab R13 and above
		S.(varglist{i}) = varglist{i+1};
	else
		warning_wrap('getargs:unknown_param', ...
		        ['Unknown parameter "' varglist{i} '"']);
	end
	i = i+2;
end
