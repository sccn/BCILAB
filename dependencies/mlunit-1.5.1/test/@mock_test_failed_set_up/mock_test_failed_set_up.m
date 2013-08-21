function m = mock_test_failed_set_up(name)
%mock_test_failed_set_up is a mock test_case with a broken set_up used for 
%the tests in test_test_case.
%
%  Example
%  =======
%         run(gui_test_runner, 'mock_test_failed_set_up(''test_method'')');

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: mock_test_failed_set_up.m 40 2006-06-11 18:24:31Z thomi $

m.dummy = 0;
w = mock_test(name);
m = class(m, 'mock_test_failed_set_up', w);


