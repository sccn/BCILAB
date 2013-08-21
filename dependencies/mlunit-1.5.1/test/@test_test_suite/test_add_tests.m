function self = test_add_tests(self)
%test_test_suite/test_add_tests tests the method test_suite/add_tests.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_suite(''test_add_tests'');');
%
%  See also test_add_tests.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_add_tests.m 47 2006-06-11 19:26:32Z thomi $

tests{1} = mock_test('test_method');
tests{2} = mock_test('test_broken_method');
self.suite = add_tests(self.suite, tests);
[self.suite, self.result] = run(self.suite, self.result);
assert(strcmp('test_result run=2 errors=1 failures=0', summary(self.result)));