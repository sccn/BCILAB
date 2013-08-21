function self = test_method(self)
%mock_test/test_method is a mock test method, that adds to the member 
%variable log the string 'test_method '.
%
%  Example
%  =======
%         test = mock_test('test_method');
%         test = run(test, self.result);
%         assert(strcmp(get_log(test), 'set_up test_method tear_down '));
%
%  See also MOCK_TEST, TEST_TEST_CASE.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_method.m 39 2006-06-11 18:15:59Z thomi $

self.log = [self.log, 'test_method '];
