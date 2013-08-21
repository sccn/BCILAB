function errors = get_error_list(self)
%test_result/get_error_list returns a cell array of tests and errors.
%
%  Example
%  =======
%  get_error_list is called for example from text_test_result/print_errors:
%         get_error_list(self)
%
%  See also TEXT_TEST_RESULT/PRINT_ERRORS.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: get_error_list.m 30 2006-06-11 15:53:00Z thomi $

errors = self.errors;