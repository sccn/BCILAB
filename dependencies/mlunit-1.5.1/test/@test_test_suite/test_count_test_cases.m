function self = test_count_test_cases(self)
%test_test_suite/test_count_test_cases tests the method
%test_suite/count_tests_cases.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_suite(''test_count_test_cases'');');
%
%  See also TEST_SUITE/COUNT_TEST_CASES.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_count_test_cases.m 47 2006-06-11 19:26:32Z thomi $

suite = test_suite;
assert(0 == count_test_cases(suite));
suite = add_test(suite, mock_test('test_method'));
assert(1 == count_test_cases(suite));
suite = add_test(suite, mock_test('test_broken_method'));
assert(2 == count_test_cases(suite));