function log = get_log(self)
%get_log returns the member variable log.
%
%  Example
%  =======
%         log = get_log(self);
%
%  See also TEST_TEST_CASE/TEST_FAILED_RESULT,
%           TEST_TEST_CASE/TEST_FAILED_SET_UP,
%           TEST_TEST_CASE/TEST_FAILED_TEAR_DOWN,
%           TEST_TEST_CASE/TEST_TEMPLATE_METHOD.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: get_log.m 39 2006-06-11 18:15:59Z thomi $

log = self.log;