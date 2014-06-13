function caller_name = determine_caller
% Determine the name of the calling user function 
%
% Out:
%   CallerName : name of the user function calling this function; ignores functions in the 
%                argument system path and helper functions (starting with hlp_)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2014-03-13

% Copyright (C) Christian Kothe, SCCN, 2014, christian@sccn.ucsd.edu
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

stack = dbstack('-completenames');
caller_name = hlp_nanocache('getcaller',10,@determine_caller_internal,stack);


% determine the name of the calling function
function caller_name = determine_caller_internal(stack)
% remove all entries from the stack that are in the arguments folder or start with hlp_
tmp = mfilename('fullpath'); tmp = strrep(tmp,[filesep 'private'],'');
arg_path = tmp(1:find(tmp==filesep,1,'last')-1);
for k=1:length(stack)
    if ~strncmp(stack(k).file,arg_path,length(arg_path)) && ~strncmp(stack(k).name,'hlp_',4)
        caller_name = stack(k).name;
        return;
    end
end
caller_name = 'commandline';
