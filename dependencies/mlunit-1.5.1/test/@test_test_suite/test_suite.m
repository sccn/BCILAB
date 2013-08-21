function self = test_suite(self)
%test_test_suite/test_suite tests the basic behaviour of test_suite/run.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_suite(''test_suite'');');
%
%  See also TEST_SUITE/RUN.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_suite.m 47 2006-06-11 19:26:32Z thomi $

self.suite = add_test(self.suite, mock_test('test_method'));
self.suite = add_test(self.suite, mock_test('test_broken_method'));
[self.suite, self.result] = run(self.suite, self.result);
assert(strcmp('test_result run=2 errors=1 failures=0', summary(self.result)));