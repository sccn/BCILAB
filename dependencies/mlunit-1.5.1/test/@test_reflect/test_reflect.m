function t = test_reflect(name)
%test_reflect tests the class reflect.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_reflect');

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_reflect.m 42 2006-06-11 18:38:25Z thomi $

t.dummy = 0;
tc = test_case(name);
t = class(t, 'test_reflect', tc);