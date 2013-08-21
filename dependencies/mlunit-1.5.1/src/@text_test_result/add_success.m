function self = add_success(self, test)
%text_test_result/add_success calls the inherited method from test_result
%and writes out an '.' (verbosity == 1) or 'OK' (verbosity == 2) to
%the stream.
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
%  $Id: add_success.m 32 2006-06-11 16:00:00Z thomi $

self.test_result = add_success(self.test_result, test);
if (self.dots)
    fprintf(self.stream, '.');
elseif (self.show_all)
    fprintf(self.stream, 'OK\n');
end;
