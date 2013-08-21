function bool = onoff2bool(onoff)
%ONOFF2BOOL        Substitutes 1 for 'on' and 0 for any other string.
%   BOOL = ONOFF2BOOL(ONOFF) returns logical 1 if ONOFF is equal to the
%   string 'on' and logical 0 otherwise.
%
%   BOOL2ONOFF and ONOFF2BOOL improve readability for Matlab settings
%   (e.g., Handle Graphics properties) that use the strings 'on' and 'off'.
%
%   Example:
%      if (onoff2bool(get(handle,'Enable'))), ...
%                         is equivalent to 
%      if (strcmp(get(handle, 'Enable'), 'on')), ...
%
%   See also BOOL2ONOFF.

bool = strcmp(onoff, 'on');
