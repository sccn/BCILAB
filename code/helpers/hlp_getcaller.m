function [name,file] = hlp_getcaller(indirection)
% Find the name & file of the calling function.
% [Name,File] = hlp_getcaller(Indirection)
%
% In:
%   Indirection : optional level of indirection (default: 1)
%
% Out: 
%   Name: MATLAB name of the calling function, if any
%   File: file name of the calling function, if any
% 
% Example:
%   % in a function, get the name of the calling function (if any)
%   mycaller = hlp_getcaller;
%
%   % in a function, get the name of the function that called the calling function
%   mycallerscaller = hlp_getcaller(2);
%
% See also:
%   hlp_iscaller
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-15

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
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

stack = dbstack;
try
    if nargout > 1
        if nargin < 1
            name = stack(3).name;
            file = stack(3).file;
        else
            name = stack(2+indirection).name;
            file = stack(2+indirection).file;
        end
    else
        if nargin < 1
            name = stack(3).name;
        else
            name = stack(2+indirection).name;
        end
    end
catch
    name = 'base';
    file = '';
end