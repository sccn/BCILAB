function self = test_get_object(self)
%test_gui_test_runner/test_get_object tests the singleton object of
%gui_test_runner.
%
%  Example
%  =======
%         run(text_test_runner, 'test_gui_test_runner(''test_get_object'')');

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_get_object.m 162 2007-01-04 11:38:53Z thomi $

self.runner = run(self.runner, '', 0);
object = get_object(gui_test_runner(1));
assert_not_equals(0, get_handle(object));
assert_equals(get_handle(object), get_handle(self.runner));