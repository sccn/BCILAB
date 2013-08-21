function print_errors(self)
%text_test_result/print_errors writes the description of all errors and 
%failures to the stream.
%
%  Example
%  =======
%  print_errors is called for example from text_test_runner/run:
%         print_errors(result);
%
%  See also TEXT_TEST_RUNNER/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: print_errors.m 32 2006-06-11 16:00:00Z thomi $

if ((self.dots) || (self.show_all))
    fprintf(self.stream, '\n');
end;
print_error_list(self, 'ERROR', get_error_list(self));
print_error_list(self, 'FAIL', get_failure_list(self));
