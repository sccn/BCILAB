function self = test_fib(name)
%test_fib constructor.
%
%  Class Info / Example
%  ====================
%  The class test_fib is the fixture for all tests of the test-driven
%  Fibonacci. The constructor shall not be called directly, but through
%  a test runner.
%
%  Run a single test method:
%         Example: run(text_test_runner, test_fib('test_null'));
%  or the whole test case (notice the single quotes around test_fib!):
%         Example: run(text_test_runner, 'test_fib');

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_fib.m 12 2006-05-26 16:11:37Z thomi $

tc = test_case(name);
self = class(struct([]), 'test_fib', tc);