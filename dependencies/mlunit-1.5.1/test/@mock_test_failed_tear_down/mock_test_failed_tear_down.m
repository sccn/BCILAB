function m = mock_test_failed_tear_down(name)
%mock_test_failed_set_up is a mock test_case with a broken tear_down used 
%for the tests in test_test_case.
%
%  Example
%  =======
%         run(gui_test_runner, 'mock_test_failed_tear_down(''test_method'')');

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: mock_test_failed_tear_down.m 40 2006-06-11 18:24:31Z thomi $

m.dummy = 0;
w = mock_test(name);
m = class(m, 'mock_test_failed_tear_down', w);


