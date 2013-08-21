function print_errors(self)
%gui_test_result/print_errors creates the list and description of all errors and 
%failures and set them to the listbox of errors and failures.
%
%  Example
%  =======
%  print_errors is called for example from gui_test_result/update:
%         print_errors(self);
%
%  See also GUI_TEST_RESULT/UPDATE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: print_errors.m 28 2006-06-11 15:11:59Z thomi $

print_error_list(self, 'ERROR', get_error_list(self), 1);
print_error_list(self, 'FAIL', get_failure_list(self));
