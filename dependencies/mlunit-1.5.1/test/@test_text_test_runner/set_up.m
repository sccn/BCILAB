function self = set_up(self)
%test_text_test_runner/set_up sets up the fixture for
%test_text_test_runner. 
%
%  Example
%  =======
%         run(gui_test_runner, 'test_text_test_runner');

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: set_up.m 147 2007-01-02 17:29:44Z thomi $

self.tmp_file = fopen('text_test_result.tmp', 'w+');
self.runner = text_test_runner(self.tmp_file, 0);
