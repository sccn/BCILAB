function self = tear_down(self)
%mock_test/tear_down is a mock tear_down, that adds to the member variable 
%log the string 'tear_down '.
%
%  Example
%  =======
%         test = mock_test('test_method');
%         test = run(test, self.result);
%         assert(strcmp(get_log(test), 'set_up test_method tear_down '));
%
%  See also MOCK_TEST, TEST_TEST_CASE.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: tear_down.m 39 2006-06-11 18:15:59Z thomi $

self.log = [self.log, 'tear_down '];