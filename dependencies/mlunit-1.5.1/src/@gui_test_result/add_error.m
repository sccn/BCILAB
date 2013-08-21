function self = add_error(self, test, error)
%gui_test_result/add_error calls the inherited method from test_result
%and the update method, which updates the different gui objects.
%
%  Example
%  =======
%  add_error is usually only called by the run method of test_case, see
%  test_case/run:
%         result = add_error(result, self, stacktrace);
%
%  See also TEST_RESULT/ADD_ERROR, TEST_CASE/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: add_error.m 28 2006-06-11 15:11:59Z thomi $

self.test_result = add_error(self.test_result, test, error);
self = update(self);