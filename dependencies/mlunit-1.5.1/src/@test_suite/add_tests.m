function self = add_tests(self, tests)
%test_suite/add_test adds a cell array of tests to the test suite.
%
%  Example
%  =======
%         suite = test_suite;
%         suite = add_tests(suite, {my_test('test_foo') ...
%             my_test('test_bar')});
%         count_test_cases(suite); % Should return 2.
%
%  See also TEST_SUITE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: add_tests.m 36 2006-06-11 16:45:23Z thomi $

if (iscell(tests))
    last = length(self.tests);
    for i = 1:length(tests)
        self.tests{last + i} = tests{i};
    end;
end;
