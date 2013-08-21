function self = test_run_failed_tests(self)
%test_text_test_runner/test_run tests the method text_test_runner/run with
%failing tests.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_text_test_runner(''test_run_failed_tests'');');
%
%  See also TEXT_TEST_RUNNER/RUN.

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_run_failed_tests.m 147 2007-01-02 17:29:44Z thomi $

% verbosity = 0
suite = test_suite;
suite = add_test(suite, test_test_case('test_template_method'));
suite = add_test(suite, test_test_case('test_broken_method'));
run(self.runner, suite);
fseek(self.tmp_file, 0, -1);

ignore_lines(self.tmp_file);
ignore_lines(self.tmp_file);

line_summary = fgetl(self.tmp_file);
pos = findstr('Ran 2 tests in ', line_summary);
if (~isempty(pos))
    assert(pos(1) == 1);
else
    message = sprintf('Test result invalid, expected <Ran 2 tests>, but was <%s>.', line_summary);
    assert(0, message);
end;

line_empty = fgetl(self.tmp_file);
assert(size(line_empty, 2) == 0);

line_ok = fgetl(self.tmp_file);
assert(strcmp('FAILED (errors=1, failures=0)', line_ok));

fclose(self.tmp_file);

function ignore_lines(fid)

ignore = 1;
while (ignore)
    line = fgetl(fid);
    if ((~ischar(line)) || ...
        (~isempty(regexp(line, '----------------------------------------------------------------------', 'once'))))
        ignore = 0;
    end;
end;
