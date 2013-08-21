function onoff = bool2onoff(bool)
%BOOL2ONOFF        Substitutes 'on' for non-zero and 'off' for zero input.
%   ONOFF = BOOL2ONOFF(BOOL) returns the string 'on' if BOOL evaluates as
%   logically true or the string 'off' otherwise. 
%
%   BOOL2ONOFF and ONOFF2BOOL improve readability for Matlab settings
%   (e.g., Handle Graphics properties) that use the strings 'on' and 'off'.
%
%   Example:
%      set(handle, 'Enable', bool2onoff(X==Y))
%          is equivalent to
%      if (X==Y), e = 'on';  else, e = 'off';  set(handle, 'Enable', e);
%
%   See also ONOFF2BOOL.

if (bool), onoff = 'on';  else  onoff = 'off';  end;
