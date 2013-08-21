function self = set_up(self)
%test_test_suite/set_up sets up the fixture for test_test_suite.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_test_suite');

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: set_up.m 47 2006-06-11 19:26:32Z thomi $

self.result = test_result;
self.suite = test_suite;
