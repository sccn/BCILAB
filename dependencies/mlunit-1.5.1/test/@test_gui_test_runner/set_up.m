function self = set_up(self)
%test_gui_test_runner/set_up sets up the fixture for test_gui_test_runner.
%
%  Example
%  =======
%         run(text_test_runner, 'test_gui_test_runner');

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: set_up.m 162 2007-01-04 11:38:53Z thomi $

self.runner = gui_test_runner;