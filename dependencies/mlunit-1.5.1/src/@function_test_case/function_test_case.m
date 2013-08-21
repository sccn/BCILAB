function self = function_test_case(test_function, set_up_function, tear_down_function)
%function_test_case contructor.
%  The constructor creates an object of the class function_test_case.
%  test_function must be a Matlab function handle to function with the
%  signature: 
%           function test_function
%
%  Class Info / Example
%  ====================
%  The class function_test_case is a wrapper for single-function tests.
%  This might be functions, which still exists, but the class provides also
%  a way for simple one-file tests (which is not possible through the
%  standard way, as deriving a test from test_case always require at least
%  two files and a directory).
%
%  The prototype of such a one-file test case looks as follows:
%           function test = test_example
%
%           test = load_tests_from_mfile(test_loader);
%
%             function set_up
%             end
%    
%             function tear_down
%             end
%
%             function test_method
%                 assert(0 == sin(0));
%             end
%  The constructor of function_test_case is called implicitly with 
%  load_test_from_mfile, which searches for the methods set_up and 
%  tear_down as well as all methods starting with "test_". For each test
%  method a separate function_test_case is created.
%           Example: function_test_case(@test_method, @set_up, @tear_down);
%
%  See also TEST_CASE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: function_test_case.m 21 2006-05-26 23:10:20Z thomi $

if (nargin == 1)
    set_up_function = 0;
    tear_down_function = 0;
elseif (nargin == 2)
    tear_down_function = 0;
end;

self.test_function = test_function;
self.set_up_function = set_up_function;
self.tear_down_function = tear_down_function;
t = test_case('', 'function_test_case');
self = class(self, 'function_test_case', t);
