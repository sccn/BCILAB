function self = add_failure(self, test, failure)
%test_result/add_failure adds a failure for the test to the list of 
%failures.
%
%  Example
%  =======
%  add_failure is usually only called by the run method of test_case, see
%  test_case/run:
%         result = add_failure(result, self, errmsg);
%
%  See also TEST_CASE/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: add_failure.m 30 2006-06-11 15:53:00Z thomi $

last = size(self.failures, 1);
self.failures{last + 1, 1} = test;
self.failures{last + 1, 2} = failure;