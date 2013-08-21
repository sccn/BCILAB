function w = mock_test(name)
%mock_test is a mock test_case used for the tests in test_test_case.
%
%  Example
%  =======
%         run(gui_test_runner, 'mock_test(''test_method'')');
%
%  See also TEST_TEST_CASE.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: mock_test.m 39 2006-06-11 18:15:59Z thomi $

w.log = '';
t = test_case(name, 'mock_test');
w = class(w, 'mock_test', t);


