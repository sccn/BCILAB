function self = run_test(self)
%function_test_case/run_test calls the test_function by the function 
%handle.
%
%  Example
%  =======
%  Usually run_test (as every test method) is not called directly, but
%  through the method run. Example:
%         test = function_test_case(@() assert(0 == sin(0)));
%         [test, result] = run(test); 
%         summary(result)
%
%  See also TEST_CASE, FUNCTION_TEST_CASE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: run_test.m 33 2006-06-11 16:02:51Z thomi $

self.test_function();
