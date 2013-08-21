function self = set_up(self) %#ok
%mock_test_failed_set_up/set_up is a mock set_up method, that is broken as
%it only class error(' ').
%
%  Example
%  =======
%         test = mock_test_failed_set_up('test_method');
%         try
%             test = run(test, default_test_result(self));
%         catch
%             assert(0);
%         end;
%         assert(strcmp('', get_log(test)));
%
%  See also TEST_TEST_CASE/TEST_FAILED_SET_UP.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: set_up.m 41 2006-06-11 18:31:37Z thomi $

error(' ');