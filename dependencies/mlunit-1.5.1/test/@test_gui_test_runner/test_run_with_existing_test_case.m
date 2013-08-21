function self = test_run_with_existing_test_case(self)
%test_gui_test_runner/test_run_with_existing_test_case tests the behaviour
%of the run method for an existing test_case.
%
%
%  Example
%  =======
%         run(text_test_runner, 'test_gui_test_runner(''test_run_with_existing_test_case'')');

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_run_with_existing_test_case.m 162 2007-01-04 11:38:53Z thomi $

run(self.runner, 'mock_test');
