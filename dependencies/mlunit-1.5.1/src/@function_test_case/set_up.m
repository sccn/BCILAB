function self = set_up(self)
%function_test_case/set_up calls the set_up_function by the function 
%handle to set up the fixture and is called everytime before 
%a test is executed.
%
%  Example
%  =======
%  set_up is not called directly, but through the method run. Example:
%         test = ... % e.g. created through load_tests_from_mfile
%         [test, result] = run(test); 
%         summary(result)
%
%  See also TEST_CASE, FUNCTION_TEST_CASE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: set_up.m 33 2006-06-11 16:02:51Z thomi $

if (strcmp(class(self.set_up_function), 'function_handle'))
    self.set_up_function();
end;