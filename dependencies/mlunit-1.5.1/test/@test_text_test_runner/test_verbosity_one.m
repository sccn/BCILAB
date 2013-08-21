function self = test_verbosity_one(self)
%test_text_test_runner/test_verbosity_one tests the method
%text_test_runner/run with verbosity = 1.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_text_test_runner(''test_verbosity_one'');');
%
%  See also TEXT_TEST_RUNNER/RUN.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_verbosity_one.m 156 2007-01-03 20:04:17Z thomi $

self.runner = text_test_runner(self.tmp_file, 1);
run(self.runner, test_test_case('test_template_method'));
fseek(self.tmp_file, 0, -1);
assert_equals('.', fgetl(self.tmp_file));
