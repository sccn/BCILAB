function self = set_up(self)
%mock_test/set_up is a mock set_up, that sets the member variable log to
%'set_up '.
%
%  Example
%  =======
%         test = mock_test('test_method');
%         test = run(test, self.result);
%         assert(strcmp(get_log(test), 'set_up test_method tear_down '));
%
%  See also MOCK_TEST, TEST_TEST_CASE.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: set_up.m 39 2006-06-11 18:15:59Z thomi $


self.log = 'set_up ';