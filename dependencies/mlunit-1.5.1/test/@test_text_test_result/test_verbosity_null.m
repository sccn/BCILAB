function self = test_verbosity_null(self)
%test_text_test_result/test_verbosity_null tests the behaviour of
%text_test_result for verbosity = 0.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_text_test_result(''test_verbosity_null'');');
%
%  See also TEXT_TEST_RESULT.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_verbosity_null.m 157 2007-01-03 20:11:10Z thomi $

result = text_test_result(self.tmp_file, 0);
set_result(self, result);
fseek(self.tmp_file, 0, -1);

line = fgetl(self.tmp_file);
assert(-1 == line);
