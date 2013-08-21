function self = test_result
%test_result contructor.
%  The constructor creates an object of the class test_result.
%
%  Class Info
%  ==========
%  The class test_result collects test results of executed tests. As in
%  the other testing frameworks of the xUnit family the framework
%  differs between failure and error. A failure is raised by an assertion,
%  that means by the method assert, while an error is raised by the Matlab 
%  environment, for example through a syntax error.
%
%  Example
%  =======
%         result = test_result;
%  
%  See also ASSERT, ASSERT_EQUALS, ASSERT_NOT_EQUALS, TEXT_TEST_RESULT.

self.tests_run = 0;
self.errors = {};
self.failures = {};
self.should_stop = 0;
self = class(self, 'test_result');