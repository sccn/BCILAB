function self = add_success(self, test)
%gui_test_result/add_success calls the inherited method from test_result
%and the update method, which updates the different gui objects.
%
%  Example
%  =======
%  add_success is usually only called by the run method of test_case, see
%  test_case/run:
%         result = add_success(result, self);
%
%  See also TEST_RESULT/ADD_SUCCESS, TEST_CASE/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: add_success.m 28 2006-06-11 15:11:59Z thomi $

self.test_result = add_success(self.test_result, test);
self = update(self);