function S = strrep_multi(S, varargin)
% ModifiedString = STRREP_MULTI(OriginalString, Substitutions...) 
% Performs multiple successive strrep-like replacements on a given string
%
% In:
%   OriginalString : The original string in which substitutions shall be performed.
%
%   Substitutions... : Arguments of the form 'oldstring1','newstring1','oldstring2','newstring2', ...
%                      where each pair of arguments defines a strrep-like substitution; the
%                      substitutions will be applied successively to the original string, same as if
%                      performing multiple strrep calls in succession.
%
% Out:
%   ModifiedString : The resulting string with substitutions performed.
%
% Examples:
%   % substitute some patterns in mystr
%   mystr = strrep_multi(mystr,'%path',mypath,'%subj',mysubj,'%session',mysession);
%
% See also: 
%   strrep
% 
%                           Christian Kothe, Syntrogi
%                           2016-04-28

% Copyright (C) Christian Kothe, Syntrogi, 2016, christian.kothe@qusp.io
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

if iscell(S)
    error('This version will only operate on a single original string.'); end

if mod(length(varargin),2) ~= 0
    error('The arguments following the original string must be even (oldstr,newstr,oldstr,newstr, ...).'); end
    
if ~iscellstr(varargin)
    error('THe arguments following the original string must all be strings.'); end

for k=1:2:length(varargin)
    S = strrep(S, varargin{k}, varargin{k+1}); end
