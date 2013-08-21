function self = test_run(self)
%test_test_case/test_result tests the method test_case/run and the
%method test_result/get_tests_run.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_case(''test_run'');');
%
%  See also TEST_CASE/RUN, TEST_RESULT/GET_TESTS_RUN.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_run.m 44 2006-06-11 18:54:09Z thomi $

test = mock_test('test_method');
[test, result] = run(test); %#ok
assert_equals(1, get_tests_run(result));
