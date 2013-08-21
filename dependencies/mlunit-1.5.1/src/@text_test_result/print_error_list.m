function print_error_list(self, prefix, errors)
%text_test_result/print_error_list is a helper function for 
%text_test_result/print_errors.
%  It iterates through all errors in errors and writes them to the stream. 
%  prefix is a string, which is written before each error or failure, to 
%  differ between them.
%
%  Example
%  =======
%  print_error_list is called twice in text_test_result/print_errors, e.g.
%  for the list of failures:
%         print_error_list(self, 'FAILURE', get_failure_list(self));
%
%  See also text_TEST_RESULT/PRINT_ERRORS.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: print_error_list.m 83 2006-10-10 19:06:10Z thomi $

for i = 1:size(errors, 1)
    for j = 1:70
        fprintf(self.stream, '=');
    end;
    fprintf(self.stream, '\n%s: %s\n', ...
        prefix, ... 
        get_description(self, errors{i, 1}));
    for j = 1:70
        fprintf(self.stream, '-');
    end;
    fprintf(self.stream, '\n%s\n', errors{i, 2});
end;
