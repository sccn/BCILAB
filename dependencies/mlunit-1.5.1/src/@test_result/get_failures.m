function failures = get_failures(self)
%test_result/get_errors returns the number of failures.
%
%  Example
%  =======
%  get_error_list is called for example from text_test_result/run:
%         get_errors(self)
%
%  See also TEXT_TEST_RESULT/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: get_failures.m 30 2006-06-11 15:53:00Z thomi $

failures = size(self.failures, 1);