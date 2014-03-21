function s = hlp_handleerror(e,level)
% Displays a formatted error message for some error object, including a full stack trace.
% hlp_handleerror(Error, Indent)
%
% In:
%   Error   : error object, as received from lasterror or via a catch clause
%             (if omitted, lasterror is used)
%   Indent  : optional indentation level, in characters
%
% Out:
%   Formatted : optionally the formatted error report
%
% Example:
%   % display an error message including stack trace in case of an error and continue without 
%   % terminating
%   try
%     ...
%   catch e
%     hlp_handleerror(e);
%   end
%
% See also:
%   lasterror, MException, catch
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-22

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

try
    if ~exist('e','var')
        e = lasterror; end %#ok<LERR>
    
    % compute the appropriate indentation level
    if ~exist('level','var')
        level = '';
    else
        level = repmat(' ',1,level);
    end
    
    if nargout == 0    
        % display the message
        for message = hlp_split(e.message,[10 13])
            fprintf('%s %s\n',level,message{1}); end
        fprintf('%s occurred in: \n',level);
        for st = e.stack'
            if ~isdeployed
                try
                    fprintf('%s   <a href="matlab:opentoline(''%s'',%i)">%s</a>: %i\n',level,st.file,st.line,st.name,st.line);
                catch
                    fprintf('%s   <a href="matlab:edit %s">%s</a>: %i\n',level,st.file,st.name,st.line);
                end
            else
                % links are not supported in deployed mode
                fprintf('%s   %s: %i\n',level,st.file,st.line);
            end
        end
    else
        % print message into a string
        s = [];
        for message = hlp_split(e.message,[10 13])
            s = [s sprintf('%s %s\n',level,message{1})]; end
        s = [s sprintf('%s occurred in: \n',level)];
        for st = e.stack'
            s = [s sprintf('%s   %s: %i\n',level,st.file,st.line)]; end
    end
catch err
    disp('An error occurred, but the traceback could not be displayed due to another error...');
end
