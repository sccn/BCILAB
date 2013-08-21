function t = test_test_case(name)
%test_test_case tests the class test_case.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_case');

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_test_case.m 44 2006-06-11 18:54:09Z thomi $

t.result = 0;
tc = test_case(name);
t = class(t, 'test_test_case', tc);