function self = test_get_methods(self)
%test_reflext/test_get_methods tests the method test_reflect/get_methods.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_reflect(''test_get_methods'')');
%
%  See also REFLECT/GET_METHODS.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_get_methods.m 47 2006-06-11 19:26:32Z thomi $

r = reflect('test_case');
m = get_methods(r);
assert(size(m, 1) > 0);
assert(~sum(strcmp(m, 'test_case')));
assert(sum(strcmp(m, 'run')) == 1);
assert(sum(strcmp(m, 'set_up')) == 1);
assert(sum(strcmp(m, 'tear_down')) == 1);
