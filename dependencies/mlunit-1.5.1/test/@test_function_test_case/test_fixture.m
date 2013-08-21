function self = test_fixture(self)
%test_function_test_case/test_fixture tests the fixture of a
%function_test_case with all test methods within on .m-file.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_function_test_case(''test_fixture'')');
%
%  See also TEST_FUNCTION_TEST_CASE.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_fixture.m 41 2006-06-11 18:31:37Z thomi $

test = function_test_case(@() assert(1), 0, 0);
[test, result] = run(test);  %#ok
assert_equals(1, get_tests_run(result));

global x; %#ok
suite = load_tests_from_mfile(test_loader);
result = test_result;
[suite, result] = run(suite, result);  %#ok
assert_equals(2, get_tests_run(result));

function set_up  %#ok

global x;
assert_equals(0, x);
x = [3 4];

function tear_down %#ok

global x;
x = 0;

function test_norm %#ok

global x;
assert_equals(5, norm(x));

function test_normest %#ok

global x;
assert_equals(5, norm(x));
