function t = test_text_test_result(name)
%test_text_test_result tests the class text_test_result.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_text_test_result');
%
%  See also TEXT_TEST_RESULT.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_text_test_result.m 157 2007-01-03 20:11:10Z thomi $

t.runner = [];
t.tmp_file = '';
tc = test_case(name, 'test_text_test_result');
t = class(t, 'test_text_test_result', tc);