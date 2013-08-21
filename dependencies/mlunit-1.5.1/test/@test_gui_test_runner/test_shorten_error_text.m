function self = test_shorten_error_text(self)
%test_gui_test_runner/test_shorten_error_text tests the method
%shorten_error_text, which shortens the directory paths of an error
%message.
%
%  Example
%  =======
%         run(text_test_runner, 'test_gui_test_runner(''test_shorten_error_text'')');

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_shorten_error_text.m 254 2007-01-27 21:23:39Z thomi $

self.runner = set_shorten(self.runner, 1);

% Empty Text
shorten_error_text(self.runner, '');

% Just one line
shorten_error_text(self.runner, 'AssertionError: No text.');

% Error text with more lines
error_text = sprintf('%s\n%s\n%s\n%s\n%s', ...
    'Traceback (most recent call first)\n', ...
    '  In s:\project\tests\mlunit\test\test_misc.m at line 18', ...
    '  In s:\project\src\@function_test_case\run_test.m at line 21', ...
    '  In s:\project\src\@test_case\run.m at line 38', ...
    'AssertionError: No text.');

shortened_text = shorten_error_text(self.runner, error_text);

error_lines = strread(error_text, '%s', 'delimiter', '\n');
shortened_lines = strread(shortened_text, '%s', 'delimiter', '\n');

assert_equals(5, length(shortened_lines));
assert_not_equals(sprintf('\n'), shortened_text(1));
assert_equals(char(error_lines{1}), char(shortened_lines{1}));
assert_equals(char(error_lines{5}), char(shortened_lines{5}));

assert_equals('In s:/../test/test_misc.m at line 18', ...
    char(shortened_lines{2}));
assert_equals('In s:/../@function_test_case/run_test.m at line 21', ...
    char(shortened_lines{3}));
assert_equals('In s:/../@test_case/run.m at line 38', ...
    char(shortened_lines{4}));
