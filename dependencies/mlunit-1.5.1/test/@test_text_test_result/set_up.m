function self = set_up(self)
%test_text_test_result/set_up sets up the fixture for
%test_text_test_result. 
%
%  Example
%  =======
%         run(gui_test_runner, 'test_text_test_result');

%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: set_up.m 157 2007-01-03 20:11:10Z thomi $

self.tmp_file = fopen('text_test_result.tmp', 'w+');

