function self = test_function_test_case(name)
%test_function_test_case tests the class function_test_case.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_function_test_case');

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_function_test_case.m 41 2006-06-11 18:31:37Z thomi $

self.dummy = 0;
test = test_case(name);
self = class(self, 'test_function_test_case', test);