function finisher = dp(message,varargin)
% stack-aware debug print
% Handle = dp(Message,Arguments...)
%
% This function allows to debug the progression of complex recursive functions. It prints a
% properly indented message immediately, and a complementing one when the function exits (or
% alternatively, when the returned handle is deleted). The two messages are of the form:
%
%   enter <functionname>: message
%   leave <functionname>: message
%
% To enable debug printout, you need to set the global variable tracking.debug.stack to true.
%
% In:
%   Message : optional message to display for the given function call (as in sprintf)
%
%   Arguments ; optional arguments to substitute inside Messsage (as in sprintf)
%
% Out:
%   Handle : Optional handle whose lifetime determines when the leave message is displayed
%
% Notes: 
%  If dp is used multiple times in the same function scope, the leave message will be displayed 
%  immediately before the next dp's enter message.
% 
% Dependencies: onCleanup()
%
% Examples:
%
%   function myfunction(...)
%   ...
%   dp('blah!');  % --> prints "enter myfunction: blah!" here
%   ...
%   < potentially recursive calls... >
%   ...
%                 % --> prints "leave myfunction: blah!" here
%
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-02-11

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

global tracking;
if isfield(tracking,'debug') && isfield(tracking.debug,'stack') && tracking.debug.stack
    % process inputs
    if ~exist('message','var')
        message = ''; end
    if ~isempty(message)
        message = sprintf([': ' message],varargin{:}); end

    % determine indentation
    stack = dbstack;
    stackdepth = length(stack);
    indent = repmat(' ',1,stackdepth*2);

    % issue leave message for any previous dp
    if nargout == 0
        assignin('caller','dp_finalizer__',[]); end

    % determine caller
    caller = stack(2).name;

    % display enter message
    disp([indent 'enter ' caller message]);

    % create leave message printer
    finisher = onCleanup(@()disp([indent 'leave ' caller message]));

    % and associate it with the calling function, if necessary
    if nargout == 0
        assignin('caller','dp_finalizer__',finisher); end
end