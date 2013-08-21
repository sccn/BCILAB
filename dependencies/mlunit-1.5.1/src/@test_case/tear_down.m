function self = tear_down(self)
%test_case/tear_down called everytime after a test is executed for 
%cleaning up the fixture.
%
%  Example
%  =======
%  tear_down is not called directly, but through the method run. Example:
%         test = ... % e.g. created through my_test('test_foo')
%         [test, result] = run(test); 
%         summary(result)
%
%  See also TEST_CASE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: tear_down.m 33 2006-06-11 16:02:51Z thomi $