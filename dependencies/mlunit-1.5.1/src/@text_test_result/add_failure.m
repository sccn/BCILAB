function self = add_failure(self, test, failure)
%text_test_result/add_failure calls the inherited method from test_result
%and writes out an 'F' (verbosity == 1) or 'FAILURE' (verbosity == 2) to
%the stream.
%
%  Example
%  =======
%  add_failure is usually only called by the run method of test_case, see
%  test_case/run:
%         result = add_failure(result, self, errmsg);
%
%  See also TEST_RESULT/ADD_FAILURE, TEST_CASE/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: add_failure.m 32 2006-06-11 16:00:00Z thomi $

self.test_result = add_failure(self.test_result, test, failure);
if (self.dots)
    fprintf(self.stream, 'F');
elseif (self.show_all)
    fprintf(self.stream, 'FAIL\n');
end;
