function self = test_should_stop(self)
%test_test_suite/test_should_stop tests the method test_result/set_should_stop.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_suite(''test_should_stop'');');
%
%  See also TEST_RESULT/SET_SHOULD_STOP.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_should_stop.m 47 2006-06-11 19:26:32Z thomi $

self.result = set_should_stop(self.result);
self.suite = add_test(self.suite, mock_test('test_method'));
[self.suite, self.result] = run(self.suite, self.result);
assert(strcmp('test_result run=0 errors=0 failures=0', summary(self.result)));

