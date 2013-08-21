function val = log_no_w( x );
% function val = log_no_w( x );
%
% val = log(x)
%
% log_no_w does not warn anything.

% to know the msgid
% log(0)
% [msg, id] = lastwarn

warning('off', 'MATLAB:log:logOfZero');
val = log(x);
warning('on', 'MATLAB:log:logOfZero');


% Local Variables: ***
% mode: matlab ***
% End: ***
