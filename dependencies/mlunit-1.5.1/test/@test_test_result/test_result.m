function self = test_result(self)
%test_test_result/test_result tests the results of test_result.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_result(''test_result'');');
%
%  See also TEST_RESULT.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_result.m 269 2007-04-02 19:54:39Z thomi $

self.result = start_test(self.result, mock_test('test_method'));
self.result = add_success(self.result, mock_test('test_method'));
assert_equals(1, was_successful(self.result));
assert(strcmp('test_result run=1 errors=0 failures=0', summary(self.result)));
self.result = add_error(self.result, mock_test('test_method'), 'foo error');
assert(strcmp('test_result run=1 errors=1 failures=0', summary(self.result)));
self.result = add_failure(self.result, mock_test('test_method'), 'foo failure');
assert(strcmp('test_result run=1 errors=1 failures=1', summary(self.result)));
self.result = stop_test(self.result, mock_test('test_method'));
self.result = set_should_stop(self.result);
assert_equals(1, get_should_stop(self.result));