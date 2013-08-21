function self = add_error(self, test, error)
%text_test_result/add_error calls the inherited method from test_result
%and writes out an 'E' (verbosity == 1) or 'ERROR' (verbosity == 2) to
%the stream.
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
%  $Id: add_error.m 32 2006-06-11 16:00:00Z thomi $

self.test_result = add_error(self.test_result, test, error);
if (self.dots)
    fprintf(self.stream, 'E');
elseif (self.show_all)
    fprintf(self.stream, 'ERROR\n');
end;