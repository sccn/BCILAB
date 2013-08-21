function self = start_test(self, test)
%gui_test_result/start_test calls the inherited method from test_result.
%
%  Example
%  =======
%  start_test is usually called by test_case/run to signal the start of the
%  test execution to the test result:
%         result = start_test(result, self);
%
%  See also TEST_RESULT/START_TEST, TEST_CASE/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: start_test.m 29 2006-06-11 15:29:58Z thomi $

self.test_result = start_test(self.test_result, test);
