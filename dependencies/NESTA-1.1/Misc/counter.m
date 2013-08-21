function out = counter( f, x )
% y = counter( f, x )
%   returns y = f(x)
%
% n = counter()
%   returns the number of calls n made to the counter and zeros-out
%   the iteration counter
%
% This uses a persistent variable, but also
%   a global variable called COUNT_
%   that may be accessed from anywhere.
%  Changing the value of COUNT_ is harmless, since
%   it will be overridden by the persistent variable
%
% Stephen Becker, srbecker@caltech.edu, March 2009
% NESTA Version 1.1


persistent n
global COUNT_
if isempty(n), n = 0; end
if isempty(COUNT_), COUNT_ = 0; end

if nargin < 1
    out = n;
    n = 0;
    return
end

out = f(x);
n = n + 1;
COUNT_ = n;
