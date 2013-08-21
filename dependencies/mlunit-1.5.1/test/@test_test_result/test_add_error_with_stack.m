function self = test_add_error_with_stack(self)
%test_test_case/test_result_list tests the method
%test_result/add_error_with_stack.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_result(''test_add_error_with_stack'');');
%
%  See also TEST_RESULT/ADD_ERROR_WITH_STACK.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_add_error_with_stack.m 269 2007-04-02 19:54:39Z thomi $

test = mock_test('test_unbalanced_parentheses');
[test, result] = run(test, self.result);
error_list = get_error_list(result);
error_lines = strread(char(error_list(2)), '%s', 'delimiter', '\n');
assert_equals('Traceback (most recent call first): ', char(error_lines(1)));
assert(~isempty(findstr('test_unbalanced_parentheses.m at line 14', char(error_lines(2)))));
assert(~isempty(findstr('run.m at line 38', char(error_lines(3)))));
assert_equals('Error: Unbalanced or misused parentheses or brackets.', char(error_lines(end - 1)));
