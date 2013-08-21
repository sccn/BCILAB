function self = test_null(self)
%test_null checks, whether the return value of fib(0) is 0.
%
%  Example
%  =======
%  Use a test runner to run the test method:
%         Example: run(text_test_runner, test_fib('test_null'));

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_null.m 14 2006-05-26 16:15:53Z thomi $

assert_equals(0, fib(0));