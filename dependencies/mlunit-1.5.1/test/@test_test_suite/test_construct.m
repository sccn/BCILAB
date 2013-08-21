function self = test_construct(self)
%test_test_suite/test_construct tests the constructor of test_suite.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_suite(''test_construct'');');
%
%  See also TEST_SUITE.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_construct.m 47 2006-06-11 19:26:32Z thomi $

tests{1} = mock_test('test_method');
tests{2} = mock_test('test_broken_method');
suite = test_suite(tests);
[suite, self.result] = run(suite, self.result); %#ok
assert_equals(2, get_tests_run(self.result));

suite = test_suite(tests);
[suite, self.result] = run(suite, self.result); %#ok
assert_equals(4, get_tests_run(self.result));
