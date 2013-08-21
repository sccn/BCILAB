function result = default_test_result(self) %#ok
%test_case/default_test_result returns s default test_result object.
%
%  Example
%  =======
%  Usually default_test_result is used by the method test_case/run to
%  obtain a default test result. If the results of more than tests should
%  be collected within the same test result, default_test_result could be
%  called before the execution of the tests. Example:
%         test1 = my_test('test_foo1');
%         test2 = my_test('test_foo2');
%         result = default_test_result(test1);
%         [test1, result] = run(test1, result)
%         [test2, result] = run(test2, result)
%         summary(result)
%
%  See also TEST_CASE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: default_test_result.m 33 2006-06-11 16:02:51Z thomi $

result = test_result;