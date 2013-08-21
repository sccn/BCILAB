function self = test_method_exists(self)
%test_reflext/test_method_exists tests the method
%test_reflect/method_exists. 
%
%  Example
%  =======
%         run(gui_test_runner, 'test_reflect(''test_method_exists'')');
%
%  See also REFLECT/METHOD_EXISTS.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_method_exists.m 47 2006-06-11 19:26:32Z thomi $

r = reflect('test_suite');
assert(method_exists(r, 'run'));
assert(~method_exists(r, 'foo'));

r = reflect('mock_test');
assert(~method_exists(r, 'foo'));
