function self = test_suite(tests)
%test_suite contructor.
%  The constructor creates an object of the class test_suite including
%  the passed given. tests is an optional parameter with a cell array of 
%  objects inherited from test_case or test_suite.
%
%  Class Info / Example
%  ====================
%  The class test_suite is a composite class to run multiple tests. A
%  test suite is created as follows:
%  Example:
%         suite = test_suite;
%         suite = add_test(suite, my_test('test_foo'));
%         suite = add_test(suite, my_test('test_bar'));
%  or
%         loader = test_loader;
%         suite = test_suite(load_tests_from_test_case(loader, 'my_test'));
%
%  Running a test suite is done the same way as a single test. Example:
%         result = test_result;
%         [suite, result] = run(suite, result);
%         summary(result)
%
%  See also TEST_CASE, TEST_LOADER, TEST_RESULT.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_suite.m 61 2006-09-21 19:11:35Z thomi $

if (nargin == 0)
    tests = [];
end;

if (iscell(tests))
    self.tests = {};
    for i = 1:length(tests)
        self.tests{i} = tests{i};
    end;
else
    self.tests = {};
end;
self.name = '';
self = class(self, 'test_suite');