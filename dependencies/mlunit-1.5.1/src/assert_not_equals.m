function assert_not_equals(expected, actual, msg)
%assert_not_equals checks whether expected ~= actual and raises an error 
%if not. 
%
%  Info / Example
%  ==============
%  A common call to assert_not_equals looks like this:
%         Example: assert_not_equals(a, b);
%  The assertion will fail, if a is equal to b. 
%
%  In addition, a message can be specified:
%         Example: assert_equals(a, b, 'a is equal to b.');
%  The message is only used, if the assertion fails.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: assert_not_equals.m 269 2007-04-02 19:54:39Z thomi $

if (nargin == 2)
    msg = '';
end;
if (isempty(msg))
    if (isnumeric(expected) || islogical(expected))
        exp_str = num2str(expected);
    else
        exp_str = expected;
    end;
    msg_out = sprintf('Expected not equal to <%s>.', exp_str);
else
    msg_out = msg;
end;
if (isequal(actual, expected))
    fail(msg_out, 1);
end;