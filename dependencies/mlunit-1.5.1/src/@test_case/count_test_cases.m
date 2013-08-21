function count = count_test_cases(self) %#ok
%test_case/count_test_cases returns the number of test cases executed 
%by run. 
%  The default implementation of test_case returns always 1, because
%  the test_case object consists only of one test method (whereas it is
%  possible to define more than one test method within the test_case class).
%
%  Example
%  =======
%         test = my_test('test_foo');
%         count_test_cases(test);     % Returns 1
%
%  See also TEST_CASE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: count_test_cases.m 83 2006-10-10 19:06:10Z thomi $

count = 1;