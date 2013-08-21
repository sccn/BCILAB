function varargout = env_handleerror(varargin)
% Displays a formatted error message for some error object, including a full stack trace.
% env_handleerror(Error, Indent)
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
%     env_handleerror(e);
%   end
%
% See also:
%   lasterror, MException, catch
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-22


[varargout{1:nargout}] = hlp_handleerror(varargin{:});
