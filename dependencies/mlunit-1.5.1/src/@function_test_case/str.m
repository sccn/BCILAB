function s = str(self)
%function_test_case/str return a string with the method and class name of 
%the test.
%
%  Example
%  =======
%  If a test method is defined as follows
%           function test_method
%               assert(0 == sin(0));
%           end
%  and a function_test_case is created for this example
%           test = function_test_case(@test_method);
%  str will return:
%           test_method(function_test_case)
%
%  See also FUNCTION_TEST_CASE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: str.m 153 2007-01-03 19:56:45Z thomi $

s = [func2str(self.test_function), '(', class(self), ')'];