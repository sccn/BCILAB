function warn_once(varargin)
% Emit a specific warning only once (per MATLAB session).
% warn_once(MessageId,Message,Arguments...)
%
% This function behaves like the built-in warning, except that a message with the given content is
% only emitted once per MATLAB session.
%
% In:
%   MessageId: optional message ID (see built-in warning)
%
%   Message : the message content, as in sprintf
%             special case: if this is 'clear', the memory of displayed warnings will be discarded
%
%   Arguments... : list of arguments that may be substituted into Message (like in sprintf)
%
% Examples:
%   % display a warning, but just once
%   warn_once('Error opening file %s',filename);
%
% See also:
%   warning
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-02-13

persistent displayed_warnings;

if strcmp(varargin{1},'clear')
    % clear the stored messages
    displayed_warnings = [];
    return;
end
    
% determine the message content
if length(varargin) > 1 && any(varargin{1}==':') && ~any(varargin{1}==' ') && ischar(varargin{2})
    message_content = [varargin{1} sprintf(varargin{2:end})];
else
    message_content = sprintf(varargin{1:end});
end

% generate a hash of it
str = java.lang.String(message_content);
message_id = sprintf('x%.0f',str.hashCode()+2^31);

% and check if it had been displayed before
if ~isfield(displayed_warnings,message_id)
    % emit the warning
    warning(varargin{:});
    % remember to not display the warning again
    displayed_warnings.(message_id) = true;
end

