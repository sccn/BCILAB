function errmsg = printf(format,varargin)
%PRINTF            Macro to display formatted text.
%   PRINTF(FORMAT,A,...) is equivalent to DISP(SPRINTF(FORMAT,A,...)),
%   displaying formatted text to the Matlab prompt.  Thus, e.g., 
%
%            PRINTF('repeat %d.', 1);
%            DISP(SPRINTF('repeat %d.', 1));
%   and      DISP(['repeat ' num2str(1) '.']);
%
%   are all are functionally equivalent.  However, (1) SPRINTF is 
%   dramatically faster than NUM2STR and (2) SPRINTF interprets control
%   sequences such as '\n' (linefeed).
%
%   ERRMSG = PRINTF(...) optionally returns any error message produced by
%   SPRINTF, or an empty matrix if no error occurred.

[s,errmsg] = sprintf(format,varargin{:});
if (isempty(errmsg)), disp(s); end;
if (nargout == 0), clear errmsg; end;