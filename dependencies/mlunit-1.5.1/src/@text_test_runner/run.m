function result = run(self, test)
%text_test_runner/run execute the test and writes the results to a stream 
%in textual form (using TextTestResult).
%
%  Example
%  =======
%      runner = text_test_runner(1, 1);
%      run(runner, mlunit_all_tests);
%
%  See also TEXT_TEST_RUNNER, TEXT_TEST_RESULT.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: run.m 149 2007-01-02 20:18:01Z thomi $

test_name = test;

if (ischar(test))
    test = load_tests_from_test_case(test_loader, test);
end;

if (strcmp(class(test), 'double') && (isempty(test)))
    fprintf(self.stream, 'Test %s not found.\n', test_name);
    return;
end;

result = text_test_result(self.stream, self.verbosity);
t = clock;
[test, result] = run(test, result); %#ok
time = etime(clock, t);
print_errors(result);
fprintf(self.stream, '----------------------------------------------------------------------\n');
tests_run = get_tests_run(result);
if (tests_run > 1) 
    fprintf(self.stream, 'Ran %d tests in %.3fs\n', tests_run, time);
else
    fprintf(self.stream, 'Ran %d test in %.3fs\n', tests_run, time);
end;
fprintf(self.stream, '\n');
if (was_successful(result))
    fprintf(self.stream, 'OK\n');
else
    fprintf(self.stream, 'FAILED (errors=%d, failures=%d)\n', get_errors(result), get_failures(result));
end;