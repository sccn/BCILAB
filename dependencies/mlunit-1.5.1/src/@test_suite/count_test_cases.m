function count = count_test_cases(self)
%test_suite/count_test_cases returns the number of test cases executed 
%by run.
%
%  Example
%  =======
%         suite = mlunit_all_tests;
%         count_test_cases(test);     % Returns 26 with version 1.4.
%
%  See also TEST_SUITE, TEST_CASE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: count_test_cases.m 36 2006-06-11 16:45:23Z thomi $

tests = length(self.tests);
count = 0;
for i = 1:tests
    count = count + count_test_cases(self.tests{i});
end;