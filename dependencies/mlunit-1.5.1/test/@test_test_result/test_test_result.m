function t = test_test_result(name)
%test_test_result tests the class test_result.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_result');
%
%  See also TEST_RESULT.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_test_result.m 46 2006-06-11 19:20:00Z thomi $

t.result = 0;
tc = test_case(name);
t = class(t, 'test_test_result', tc);