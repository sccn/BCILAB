function self = test_failed_tear_down(self)
%test_test_case/test_failed_tear_down tests the behaviour of test_case/run,
%if the tear_down method fails.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_case(''test_failed_tear_down'');');
%
%  See also TEST_CASE/RUN.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_failed_tear_down.m 44 2006-06-11 18:54:09Z thomi $

test = mock_test_failed_tear_down('test_method');
try
    test = run(test, default_test_result(self));
catch
    assert(0);
end;
assert(strcmp('set_up test_method ', get_log(test)));
