function self = test_verbosity_null(self)
%test_text_test_runner/test_verbosity_null tests the method
%text_test_runner/run with verbosity = 0.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_text_test_runner(''test_verbosity_null'');');
%
%  See also TEXT_TEST_RUNNER/RUN.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_verbosity_null.m 156 2007-01-03 20:04:17Z thomi $

% verbosity = 0
run(self.runner, test_test_case('test_template_method'));
fseek(self.tmp_file, 0, -1);

line_sep = fgetl(self.tmp_file);
assert(strcmp('----------------------------------------------------------------------', line_sep));

line_summary = fgetl(self.tmp_file); 
pos = findstr('Ran 1 test in ', line_summary);
if (~isempty(pos))
    assert(pos(1) == 1);
else
    assert(0);
end;

assert_equals(0, size(fgetl(self.tmp_file), 2));
assert_equals('OK', fgetl(self.tmp_file));

