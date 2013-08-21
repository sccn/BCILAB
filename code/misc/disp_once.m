function got_displayed = disp_once(varargin)
% Display a specific message only once (per MATLAB session).
% GotDisplayed = disp_once(Message,Arguments...)
%
% This function displays a message like fprintf (though with implicit newline), except that a 
% message with the given content is only emitted once per MATLAB session.
%
% In:
%   Message : the message content (as in fprintf)
%             special case: if this is 'clear', the memory of displayed messages will be discarded
%
%   Arguments... : list of arguments that may be substituted into Message (as in fprintf)
%
% Out:
%   GotDisplayed : whether the message got displayed
%
% Examples:
%   % display a message, but just once
%   disp_once('Note: Trying to read file now...');
%
% See also:
%   disp
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-02-13

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

persistent displayed_messages;

if strcmp(varargin{1},'clear')
    % clear the stored messages
    displayed_messages = [];
    return;
end

message_content = sprintf(varargin{1:end});

% generate a hash of it
str = java.lang.String(message_content);
message_id = sprintf('x%.0f',str.hashCode()+2^31);

% and check if it had been displayed before
if ~isfield(displayed_messages,message_id)
    % emit the message
    disp(message_content);
    got_displayed = true;
    % remember to not display the message again
    displayed_messages.(message_id) = true;
else
    got_displayed = false;
end

