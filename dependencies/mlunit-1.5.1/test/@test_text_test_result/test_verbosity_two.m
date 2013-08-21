function self = test_verbosity_two(self)
%test_text_test_result/test_verbosity_two tests the behaviour of
%text_test_result for verbosity = 2.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_text_test_result(''test_verbosity_two'');');
%
%  See also TEXT_TEST_RESULT.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_verbosity_two.m 157 2007-01-03 20:11:10Z thomi $

result = text_test_result(self.tmp_file, 2);
result = set_result(self, result);
print_errors(result);
fseek(self.tmp_file, 0, -1);

line = fgetl(self.tmp_file);
assert(strcmp('mock_test(''test_method'') ... OK', line));
line = fgetl(self.tmp_file);
assert(strcmp('mock_test(''test_method'') ... ERROR', line));
line = fgetl(self.tmp_file);
assert(strcmp('mock_test(''test_method'') ... FAIL', line));

line = fgetl(self.tmp_file);
assert(strcmp('', line));
line = fgetl(self.tmp_file);
assert(strcmp('======================================================================', line));
line = fgetl(self.tmp_file);
assert(strcmp('ERROR: mock_test(''test_method'')', line));
line = fgetl(self.tmp_file);
assert(strcmp('----------------------------------------------------------------------', line));
line = fgetl(self.tmp_file);
assert(strcmp('foo error', line));
line_ignore = fgetl(self.tmp_file); %#ok
line_ignore = fgetl(self.tmp_file); %#ok
line = fgetl(self.tmp_file);
assert(strcmp('======================================================================', line));
line = fgetl(self.tmp_file);
assert(strcmp('FAIL: mock_test(''test_method'')', line));
line = fgetl(self.tmp_file);
assert(strcmp('----------------------------------------------------------------------', line));
line = fgetl(self.tmp_file);
assert(strcmp('foo failure', line));
