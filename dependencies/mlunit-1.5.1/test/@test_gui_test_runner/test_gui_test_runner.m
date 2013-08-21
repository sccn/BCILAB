function self = test_gui_test_runner(name)
%test_gui_test_runner tests the class gui_test_runner.
%
%  Example
%  =======
%         run(text_test_runner, 'test_gui_test_runner');

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: test_gui_test_runner.m 162 2007-01-04 11:38:53Z thomi $

self.runner = [];
tc = test_case(name);
self = class(self, 'test_gui_test_runner', tc);