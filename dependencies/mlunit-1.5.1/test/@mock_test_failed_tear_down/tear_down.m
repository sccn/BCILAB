function self = tear_down(self) %#ok
%mock_test_failed_tear_down/tear_down is a mock tear_down method, that is
%broken as it only class error(' ').
%
%  Example
%  =======
%         test = mock_test_failed_tear_down('test_method');
%         try
%             test = run(test, default_test_result(self));
%         catch
%             assert(0);
%         end;
%         assert(strcmp('set_up test_method ', get_log(test)));
%
%  See also TEST_TEST_CASE/TEST_FAILED_TEAR_DOWN.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: tear_down.m 41 2006-06-11 18:31:37Z thomi $

error(' ');