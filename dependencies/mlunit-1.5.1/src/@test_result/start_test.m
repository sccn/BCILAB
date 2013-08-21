function self = start_test(self, test) %#ok
%test_result/start_test indicates that a test will be started.
%
%  Example
%  =======
%  start_test is usually called by test_case/run to signal the start of the
%  test execution to the test result:
%         result = start_test(result, self);
%
%  See also TEST_CASE/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: start_test.m 30 2006-06-11 15:53:00Z thomi $

self.tests_run = self.tests_run + 1;