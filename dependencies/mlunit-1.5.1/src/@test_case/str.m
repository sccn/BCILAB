function s = str(self)
%test_case/str returns a string with the method and class name of the test.
%
%  Example
%  =======
%  If a test method is defined as follows
%           function test_method
%               assert(0 == sin(0));
%           end
%  belonging to a class my_test, which is instantiated for example
%           test = my_test('test_method');
%  str will return:
%           my_test('test_method')
%
%  See also TEST_CASE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: str.m 61 2006-09-21 19:11:35Z thomi $

s = [class(self), '(''', self.name, ''')'];