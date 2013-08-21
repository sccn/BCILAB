function s = summary(self)
%test_result/summary returns a string with a summary of the test result.
%
%  Example
%  =======
%         test = ... % e.g. created through my_test('test_foo')
%         [test, result] = run(test); 
%         summary(result)

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: summary.m 30 2006-06-11 15:53:00Z thomi $

s = sprintf('%s run=%d errors=%d failures=%d', class(self), self.tests_run, get_errors(self), get_failures(self));
