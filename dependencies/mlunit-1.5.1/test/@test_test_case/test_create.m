function self = test_create(self)
%test_test_case/test_create tests the constructor of test_case.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_case(''test_creates'');');
%
%  See also TEST_CASE.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_create.m 44 2006-06-11 18:54:09Z thomi $

error = 0;
try
    test_case('foo', 'mock_test');
    error = 1;
catch
end;
assert_equals(0, error);

error = 0;
try
    test_case('', 'mock_test');
catch
end;
assert_equals(0, error);

