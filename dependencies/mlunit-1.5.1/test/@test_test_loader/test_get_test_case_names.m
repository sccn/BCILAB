function self = test_get_test_case_names(self)
%test_test_loader/test_get_test_case_names tests the method
%test_loader/get_test_case_names.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_loader(''test_get_test_case_names'');');
%
%  See also TEST_LOADER/GET_TEST_CASE_NAMES.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_get_test_case_names.m 45 2006-06-11 18:58:55Z thomi $

t = test_loader;
n = get_test_case_names(t, 'mock_test');
assert(size(n, 1) > 0);
assert(sum(strncmp(n, 'test', 4)) == size(n, 1));

n = get_test_case_names(t, 'mock_test_failed_set_up');
% Number of methods is zero as inheritance is not supported by reflect
assert(sum(strcmp(n, 'test_method')) == 0);