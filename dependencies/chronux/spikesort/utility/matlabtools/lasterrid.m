function lastid = lasterrid
%LASTERRID         Message ID associated with last error.
%   LASTID = LASTERRID returns a string containing the identifier string
%   corresponding to the last error message issued by MATLAB (see HELP
%   ERROR for more information on message identifiers).
%
%   See also LASTERR.

[lastmsg, lastid] = lasterr;
