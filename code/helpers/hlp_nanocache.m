function varargout = hlp_nanocache(groupname,poolsize,func,varargin)
% Cache results of function calls.
% Results... = hlp_nanocache(Groupname,Poolsize,Function,Arguments...)
%
% The nanocache is useful when you want to cache the result of a computation for later reuse. This
% is done by comparing the given inputs with cached copies using isequal(), so you can only have a
% very small number of distinct inputs unless the inputs differ in their byte size (like data
% structures), in which case the nanocache is extremely efficient. You can also map them onto
% different groupnames (e.g., using a hash function), in which case the nanocache is also very
% efficient.
%
% In:
%   Groupname : A name that separates different groups of nanocache calls (e.g., name of the calling
%               function, or name of a key string in the Arguments). Can also be the same for every call.
%
%   Poolsize : Maximum number of entries to store per cache pool. Each group has one cache pool for
%              each Function and each length of Arguments (in bytes). The pool size limits how many 
%              arguments of the same length the nanocache can hold for a given Function and Groupname.
%              Recommended number: 10.
%
%   Function : Function handle or name that should be called.
%
%   Arguments... : Arguments to pass to the function.
%
% Out:
%   Results... : Return values of the function for the given arguments, or cached copy.
%
% Limitations:
%   You cannot use the nanocache for functions whose output depends on anything other than the
%   Arguments (like global variables, files on disk, random numbers, and so on). However, you can make
%   each of these an extra argument to the function in order to use the nanocache. 
%
%   If a function should give different results even through the arguments have the same value under
%   isequal, for instance when the type differs but the value is the same (like uint32(3) vs.
%   int32(3)), you need to make the type an extra argument to the function.
%
%   If the Arguments contain NaNs or anonymous functions the nanocache is not effective.
%
% Notes:
%   There is another function that supports more features to overcome most of the limitations of the
%   nanocache (but with extra overhead), called hlp_microcache.
%
%   To clear the nanocache, run the command "clear hlp_nanocache"
%
% Examples:
%   % if this line is executed for the first time, it is as slow as num2str(20.5)
%   str = hlp_nanocache('conversion',10,@num2str,20.5);
%
%   % if it is executed a second time, it is more than 60% faster than num2str(20.5)
%   str = hlp_nanocache('conversion',10,@num2str,20.5);
%
%                                   Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                   2013-10-16

% Copyright (C) Christian Kothe, SCCN, 2013, christian@sccn.ucsd.edu
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

    persistent cache;
    
    % get the pool id based on the call signature (combination of nargout and input size)
    varinfo = whos('varargin');
    poolid = [char('b'+nargout) sprintf('%u',varinfo.bytes)];
    
    try
        % retrieve the pool of size-equivalent objects
        pool = cache.(groupname).(char(func)).(poolid);
        % search for the key in the pool
        for k=1:length(pool.keys)
            if isequal(varargin,pool.keys{k})
                varargout = pool.values{k};
                return;
            end
        end
    catch %#ok<CTCH>
        % create new pool
        pool = struct('keys',{{}},'values',{{}});
    end
    
    % did not find the entry: compute it
    [varargout{1:nargout}] = feval(func,varargin{:});
    
    % add to pool (max size per pool is limited)
    idx = 1+mod(length(pool.keys),poolsize);
    pool.keys{idx} = varargin;
    pool.values{idx} = varargout;
    cache.(groupname).(char(func)).(poolid) = pool;
end
