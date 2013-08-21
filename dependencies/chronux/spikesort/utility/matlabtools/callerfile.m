function caller = callerfile(option)
%CALLERFILE        Name of calling M-file.
%   When called from within an M-file, CALLERFILE returns a string
%   containing the name of the preceding M-file on the call stack.  If the
%   M-file was invoked directly from the command line (or if CALLERFILE is 
%   itself invoked directly from the command line), CALLERFILE returns an
%   empty string.
%
%   P = CALLERFILE('fullpath') returns the full path and name of the
%   M-file preceding the M-file in which the call occurs, without the
%   extension.
%
%   If CALLERFILE is called with any argument other than 'fullpath', it
%   behaves as if it were called with no argument.
%
%   See also DBSTACK, INPUTNAME, MFILENAME.

stack = dbstack;     % get the call stack

if (length(stack) > 2)   % if it goes beyond the requesting m-file
    [p,caller] = fileparts(stack(3).name);
    if ((nargin > 0) && (strcmp(option, 'fullpath')))
        caller = [p caller];
    end
else
    caller = '';    
end

