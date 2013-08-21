function self = test_run_with_nonexisting_test_case(self)
%test_text_test_runner/test_run_with_nonexisting_test_case tests the
%behaviour of the run method for a nonexisting test_case.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_text_test_runner(''test_run_with_nonexisting_test_case'');');
%
%  See also TEXT_TEST_RUNNER/RUN.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_run_with_nonexisting_test_case.m 147 2007-01-02 17:29:44Z thomi $

run(self.runner, 'mlunit_nonexisting_test');
fseek(self.tmp_file, 0, -1);
assert_equals('Test mlunit_nonexisting_test not found.', fgetl(self.tmp_file));
