function self = add_error(self, test, error)
%test_result/add_error adds an error for the test to the list of errors.
%
%  Example
%  =======
%  add_error is usually only called by the run method of test_case, see
%  test_case/run:
%         result = add_error(result, self, stacktrace);
%
%  See also TEST_CASE/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: add_error.m 30 2006-06-11 15:53:00Z thomi $

newlines = sprintf('\n\n');
if (size(strfind(error, newlines)) == 0)
    error = sprintf('%s\n\n', error);
end;
last = size(self.errors, 1);
self.errors{last + 1, 1} = test;
self.errors{last + 1, 2} = error;