function texts(self)
%gui_test_result/texts creates the text area with the number of tests ran, the
%number of errors and the number of failures.
%
%  Example
%  =======
%  texts is called for example from gui_test_result/update:
%         texts(self);
%
%  See also GUI_TEST_RESULT/UPDATE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: texts.m 28 2006-06-11 15:11:59Z thomi $

set(self.text_runs, 'String', ['Runs: ', num2str(get_tests_run(self)), ...
    ' / Errors: ', num2str(get_errors(self)), ...
    ' / Failures: ', num2str(get_failures(self))]);
