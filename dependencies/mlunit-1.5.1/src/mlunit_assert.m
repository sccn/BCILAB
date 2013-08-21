function assert(expr, msg, internal_call)
%assert checks whether the expr is true or not and raises an error if not. 
%
%  Info / Example
%  ==============
%  A common call to assert looks like this:
%         Example: assert(a == b);
%  The assertion will fail, if a is not equal to b. 
%
%  In addition, a message can be specified:
%         Example: assert(a == b, 'a is not equal to b.');
%  The message is only used, if the assertion fails.
%
%  The parameter internal_call is used to trim the top entries of the stack
%  trace, e.g. for user-defined assert methods (see the code of
%  assert_equals for an example).

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: assert.m 276 2007-04-06 15:32:12Z thomi $

if (nargin == 1)
    msg = 'no message.';
    internal_call = 0;
end;
if (nargin == 2)
    internal_call = 0;
end;
if (internal_call > 1)
    internal_call = 1;
end;
if ((isempty(expr)) || (~expr))
    stack = dbstack('-completenames');
    stacktrace = '';
    for i = 2 + internal_call:size(stack, 1)
        stacktrace = sprintf('%s\n  In %s at line %d', ...
            stacktrace, ...
            stack(i).file, stack(i).line);
    end;
    stacktrace = sprintf('%s\n', stacktrace);
    error(['MLUNIT FAILURE:Traceback (most recent call first): ', ...
        stacktrace, ...
        'AssertionError: ', ...
        msg]);
end;
