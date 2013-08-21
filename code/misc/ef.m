function ef(varargin)
% Command-style shortcut for exp_fullform(varargin)
%
% See also: 
%   exp_fullform
%
% Example:
% >> ef var1 var2
%    var1 = 
%      ...
%    var2 = 
%      ...
%

% Copyright (C) Christian Kothe, SCCN, 2011, christian@sccn.ucsd.edu
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

if ~iscellstr(varargin)
    error('ef is to be used as a command, i.e. without brackets.'); end

for i=1:length(varargin)
    disp([varargin{i} ' = ']);
    disp(['  ' exp_fullform(evalin('caller',varargin{i}))]);
    disp('');
end