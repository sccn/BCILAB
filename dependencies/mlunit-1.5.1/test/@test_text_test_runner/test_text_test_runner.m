function t = test_text_test_runner(name)
%test_text_test_runner tests the class text_test_runner.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_text_test_runner');
%
%  See also TEXT_TEST_RUNNER.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_text_test_runner.m 156 2007-01-03 20:04:17Z thomi $

t.runner = [];
t.tmp_file = 0;
tc = test_case(name, 'test_text_test_runner');
t = class(t, 'test_text_test_runner', tc);