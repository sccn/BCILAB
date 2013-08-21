function self = stop_test(self, test) %#ok
%test_result/stop_test indicates that a test has been finished.
%
%  Example
%  =======
%  stop_test is usually called by test_case/run to signal the end of the
%  test execution to the test result:
%         result = stop_test(result, self);
%
%  See also TEST_CASE/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: stop_test.m 30 2006-06-11 15:53:00Z thomi $

