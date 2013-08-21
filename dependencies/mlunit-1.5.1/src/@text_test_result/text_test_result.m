function self = text_test_result(stream, verbosity)
%text_test_result contructor.
%  The constructor creates an object of the class text_test_result. The
%  output is written to the value of stream. The parameter verbosity
%  defines, how much output is written. Possible values are 0, 1 and 2.
%
%  Class Info
%  ==========
%  The class text_test_result is inherited from test_result and prints
%  formatted test results to a stream.
%
%  Example
%  =======
%  Output all results to the Matlab Command Window:
%         result = text_test_result(1, 0)
%
%  See also TEST_RESULT.

self.stream = stream;
if (verbosity == 1)
    self.dots = 1;
    self.show_all = 0;
elseif (verbosity > 1)
    self.dots = 0;
    self.show_all = 1;
else
    self.dots = 0;
    self.show_all = 0;
end;
self.verbosity = verbosity;
result = test_result;
self = class(self, 'text_test_result', result);