function self = add_error_with_stack(self, test, err)
%test_result/add_error_with_stack adds an error for the test to the list of
%errors using the struct members of lasterror (file and stack) to generate
%a more detailed error report.
%
%  Example
%  =======
%         result = add_error_with_stack(result, self, lasterror);
%
%  See also TEST_RESULT/ADD_ERROR, TEST_CASE/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: add_error_with_stack.m 249 2007-01-26 22:59:59Z thomi $

stacktrace = '';
for i = 1:size(err.stack, 1)
    stacktrace = sprintf('%s\n  In %s at line %d', ...
        stacktrace, ...
        err.stack(i).file, err.stack(i).line);
end;
[message, stacktrace] = parse_error(self, err.message, stacktrace);
stacktrace = sprintf('%s\n', stacktrace);
self = add_error(self, ...
    test, ...
    ['Traceback (most recent call first): ', ...
    stacktrace, ...,
    'Error: ', ...
    message]);