function self = test_broken_method(self) %#ok
%mock_test/test_broken_method is a mock test method, that is broken as it
%only class error(' ').
%
%  Example
%  =======
%         test = mock_test('test_broken_method');
%         [test, self.result] = run(test, self.result);
%         assert_equals('test_result run=1 errors=1 failures=0', summary(self.result));
%         assert(strcmp('set_up tear_down ', get_log(test)));
%
%  See also MOCK_TEST, TEST_TEST_CASE.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_broken_method.m 41 2006-06-11 18:31:37Z thomi $
error(' ');