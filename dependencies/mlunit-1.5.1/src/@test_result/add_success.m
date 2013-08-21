function self = add_success(self, test) %#ok
%test_result/add_success is an empty method for classes, which might
%do some action on a successful test.
%
%  Example
%  =======
%  add_success is usually only called by the run method of test_case, see
%  test_case/run:
%         result = add_success(result, self);
%
%  See also TEXT_TEST_RESULT/ADD_SUCCESS, TEST_CASE/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: add_success.m 30 2006-06-11 15:53:00Z thomi $

