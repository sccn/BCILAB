function self = add_test(self, test)
%test_suite/add_test adds a test to the test suite.
%  If test is empty, nothing is done.
%
%  Example
%  =======
%         suite = test_suite;
%         suite = add_test(suite, my_test('test_foo'));
%         suite = add_test(suite, my_test('test_bar'));
%         count_test_cases(suite); % Should return 2.
%
%  See also TEST_SUITE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: add_test.m 80 2006-10-09 19:15:48Z thomi $

if (~isempty(test))
    self.tests{length(self.tests) + 1} = test;
end;
