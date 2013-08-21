function self = test_count_test_cases(self)
%test_test_case/test_count_test_cases tests the method
%test_case/count_test_cases, whose return value has to be one.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_case(''test_count_test_cases'');');
%
%  See also TEST_CASE/COUNT_TEST_CASES.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_count_test_cases.m 44 2006-06-11 18:54:09Z thomi $

assert(1 == count_test_cases(self));