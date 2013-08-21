function self = test_value(self)
%test_value tests different values of the fibonacci function (y = fib(x)).
%
%  Example
%  =======
%  Use a test runner to run the test method:
%         Example: run(text_test_runner, test_fib('test_value'));

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_value.m 15 2006-05-26 16:17:55Z thomi $

assert_equals(1, fib(1));
assert_equals(1, fib(2));
assert_equals(2, fib(3));
assert_equals(3, fib(4));
assert_equals(5, fib(5));
assert_equals(8, fib(6));
assert_equals(13, fib(7));
assert_equals(21, fib(8));
assert_equals(34, fib(9));
assert_equals(55, fib(10));