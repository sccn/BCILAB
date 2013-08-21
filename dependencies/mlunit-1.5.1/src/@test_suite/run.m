function [self, result] = run(self, result)
%test_suite/run executes the test suite and saves the results in result.
%
%  Example
%  =======
%  Running a test suite is done the same way as a single test. Example:
%         suite = ...; % Create test_suite, e.g. with test_loader.
%         result = test_result;
%         [suite, result] = run(suite, result);
%         summary(result)
%
%  See also TEST_SUITE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: run.m 36 2006-06-11 16:45:23Z thomi $
for i = 1:length(self.tests)
    if (get_should_stop(result))
        break;
    end;

    [self.tests{i}, result] = run(self.tests{i}, result);
end;
    