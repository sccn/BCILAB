function self = text_test_runner(stream, verbosity)
%text_test_runner contructor.
%  The constructor creates an object of the class text_test_runner.
%
%  Class Info / Example
%  ====================
%  The class text_test_runner runs a test_case or test_suite and writes the 
%  results to a stream in textual form (using text_test_result). 
%  Example:
%      runner = text_test_runner(1, 1);
%      run(runner, mlunit_all_tests);
%
%  See also TEXT_TEST_RESULT.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: text_test_runner.m 34 2006-06-11 16:07:51Z thomi $

if (nargin == 0)
    stream = 1;
    verbosity = 0;
end;

self.stream = stream;
self.verbosity = verbosity;
self = class(self, 'text_test_runner');