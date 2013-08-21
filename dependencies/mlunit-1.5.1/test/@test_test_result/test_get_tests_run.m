function self = test_get_tests_run(self)
%test_test_result/test_get_tests_run tests the method
%test_result/get_tests_run.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_result(''test_get_tests_run'');');
%
%  See also TEST_RESULT/GET_TESTS_RUN.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_get_tests_run.m 46 2006-06-11 19:20:00Z thomi $

self.result = start_test(self.result, mock_test('test_method'));
self.result = add_success(self.result, mock_test('test_method'));
self.result = stop_test(self.result, mock_test('test_method'));
assert(1 == get_tests_run(self.result));
