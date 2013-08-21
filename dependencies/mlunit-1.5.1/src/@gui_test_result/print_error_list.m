function print_error_list(self, prefix, errors, reset_list)
%gui_test_result/print_error_list is a helper function for 
%gui_test_result/print_errors.
%  It iterates through all errors in errors, creates a cell array 
%  containing the error title and a cell array containing the description. 
%  The first is saved as String, the second as UserData of the listbox, 
%  which is shown when selecting an error.
%
%  Example
%  =======
%  print_error_list is called twice in gui_test_result/print_errors, e.g.
%  for the list of failures:
%         print_error_list(self, 'FAILURE', get_failure_list(self));
%
%  See also GUI_TEST_RESULT/PRINT_ERRORS.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: print_error_list.m 83 2006-10-10 19:06:10Z thomi $

if (nargin == 3)
    reset_list = 0;
end;

list = get(self.error_listbox, 'String');
data = get(self.error_listbox, 'UserData');
if ((isempty(list)) || (reset_list == 1))
    list = cell(0);
    data = cell(0);
end;

for i = 1:size(errors, 1)
    idx = length(list) + 1;
    list{idx} = sprintf('%s: %s', ...
        prefix, ...
        get_description(self, errors{i, 1}));
    data{idx} = errors{i, 2};
end;

set(self.error_listbox, 'String', list);
set(self.error_listbox, 'UserData', data);
set(self.error_listbox, 'Value', 1);
