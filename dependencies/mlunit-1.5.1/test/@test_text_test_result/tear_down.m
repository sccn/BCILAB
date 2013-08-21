function self = tear_down(self)
%test_text_test_result/tear_down tears down the fixture for
%test_text_test_result.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_text_test_result');

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: tear_down.m 157 2007-01-03 20:11:10Z thomi $

try
    fclose(self.tmp_file);
    delete('text_test_result.tmp');
catch
end;
