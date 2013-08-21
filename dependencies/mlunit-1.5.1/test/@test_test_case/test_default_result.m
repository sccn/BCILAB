function self = test_default_result(self)
%test_test_case/test_default_result tests the method
%test_case/default_test_result.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_case(''test_default_result'');');
%
%  See also TEST_CASE/DEFAULT_TEST_RESULT.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_default_result.m 44 2006-06-11 18:54:09Z thomi $

t = mock_test('test_method');
assert(isa(default_test_result(t), 'test_result'));