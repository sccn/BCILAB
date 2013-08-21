function t = test_assert(name)
%test_assert tests the methods assert, assert_equals and assert_not_equals.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_assert');

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_assert.m 54 2006-07-09 16:41:26Z thomi $

tc = test_case(name);
t = class(struct([]), 'test_assert', tc);