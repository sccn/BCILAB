function t = test_test_suite(name)
%test_test_suite tests the class test_suite.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_suite');
%
%  See also TEST_SUITE.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_test_suite.m 47 2006-06-11 19:26:32Z thomi $

t.result = 0;
t.suite = 0;
tc = test_case(name);
t = class(t, 'test_test_suite', tc);