function self = test_result(self)
%test_test_case/test_result tests the method test_case/run and the
%test_result. 
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_case(''test_result'');');
%
%  See also TEST_CASE/RUN, TEST_RESULT.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_result.m 44 2006-06-11 18:54:09Z thomi $

test = mock_test('test_method');
[test, self.result] = run(test, self.result); %#ok
assert(strcmp('test_result run=1 errors=0 failures=0', summary(self.result)));