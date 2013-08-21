function self = test_get_errors_failures(self)
%test_test_result/test_get_errors_failures tests the methods
%test_result/get_errors and test_result/get_failures.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_result(''test_get_errors_failures'');');
%
%  See also TEST_RESULT/GET_ERRORS, TEST_RESULT/GET_FAILURES.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_get_errors_failures.m 46 2006-06-11 19:20:00Z thomi $

self.result = start_test(self.result, mock_test('test_method'));
self.result = add_error(self.result, mock_test('test_method'), 'foo error');
self.result = add_failure(self.result, mock_test('test_method'), 'foo failure');
self.result = stop_test(self.result, mock_test('test_method'));
assert(1 == get_errors(self.result));
assert(1 == get_failures(self.result));
