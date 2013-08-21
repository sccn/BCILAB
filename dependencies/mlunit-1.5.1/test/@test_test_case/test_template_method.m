function self = test_template_method(self)
%test_test_case/test_result tests the method test_case/run.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_case(''test_template_method'');');
%
%  See also TEST_CASE/RUN.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_template_method.m 44 2006-06-11 18:54:09Z thomi $

test = mock_test('test_method');
test = run(test, self.result);
assert(strcmp(get_log(test), 'set_up test_method tear_down '));
