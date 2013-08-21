function tests_run = get_tests_run(self)
%test_result/tests_run returns the number of tests executed.
%
%  Example
%  =======
%  get_tests_run is called for example from text_test_runner/run:
%         tests_run = get_tests_run(result);
%
%  See also TEXT_TEST_RUNNER/RUN.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: get_tests_run.m 30 2006-06-11 15:53:00Z thomi $

tests_run = self.tests_run;