function self = set_up(self)
%test_case/set_up sets up the fixture and is called everytime before 
%a test is executed.
%
%  Example
%  =======
%  set_up is not called directly, but through the method run. Example:
%         test = ... % e.g. created through my_test('test_foo')
%         [test, result] = run(test); 
%         summary(result)
%
%  See also TEST_CASE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: set_up.m 33 2006-06-11 16:02:51Z thomi $