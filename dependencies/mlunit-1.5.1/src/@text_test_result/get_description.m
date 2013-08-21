function description = get_description(self, test) %#ok
%text_test_result/get_description returns the name of the test.
%
%  Example
%  =======
%  get_description is called by print_error_list to get the name of the
%  test, in which an error or failure occured. See
%  text_test_result/print_error_list:
%         get_description(self, errors{i, 1})
%
%  See also TEXT_TEST_RESULT/PRINT_ERROR_LIST.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: get_description.m 32 2006-06-11 16:00:00Z thomi $

description = str(test);