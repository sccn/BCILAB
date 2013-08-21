function id = fresh_id(tag)
% Get an integer that is fresh/unique for a given tag.
% 
% In:
%   Tag : arbitrary string tag (must conform to MATLAB variable naming rules). The generated id's
%         are independent for different tags and local to the MATLAB session/instance.
%
% Out:
%   Id : An integer id that is fresh (i.e. unique), for a particular tag. When a new tag is first
%        used, this function will return 1. On every further call it will return the next higher 
%        unused integer. Closing the MATLAB session or calling "clear all" will reset any counter
%        back to 1 (also note that ids are not guaranteed to be unique across machines that run in 
%        parallel).
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-11-24

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

persistent ids;

try
    % get next higher id for the given tag
    id = ids.(tag).incrementAndGet();
catch
    try
        % tag doesn't exist yet: create
        ids.(tag) = java.util.concurrent.atomic.AtomicInteger();
        id = ids.(tag).incrementAndGet();
    catch
        if ~exist('tag','var')
            error('Please specify a tag for which you would like to obtain an id.'); end
        if ~isvarname(tag)
            error('Tags must be valid MATLAB variable names.'); end
    end
end
