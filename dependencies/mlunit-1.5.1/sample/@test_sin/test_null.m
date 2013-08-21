function self = test_null(self)
%test_null checks, whether the return value of sin(0) is 0.
%
%  Example
%  =======
%  Use a test runner to run the test method:
%         Example: run(text_test_runner, test_sin('test_null'));

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_null.m 13 2006-05-26 16:14:26Z thomi $

assert_equals(0, sin(0));