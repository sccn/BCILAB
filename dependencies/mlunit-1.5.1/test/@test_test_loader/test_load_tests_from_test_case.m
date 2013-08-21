function self = test_load_tests_from_test_case(self)
%test_test_loader/test_load_tests_from_test_case tests the method
%test_loader/load_tests_from_test_case.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_loader(''test_load_tests_from_test_case'');');
%
%  See also TEST_LOADER/LOAD_TESTS_FROM_TEST_CASE.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_load_tests_from_test_case.m 253 2007-01-27 21:20:20Z thomi $

t = test_loader;
suite = load_tests_from_test_case(t, 'mock_test');
result = test_result;
[suite, result] = run(suite, result); %#ok
assert_equals(3, get_tests_run(result));
assert_equals(2, get_errors(result));
