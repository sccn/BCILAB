function self = tear_down(self)
%test_text_test_runner/tear_down tears down the fixture for
%test_text_test_runner.
%
%  Example
%  =======
%         run(gui_test_runner, 'test_text_test_runner');

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: tear_down.m 147 2007-01-02 17:29:44Z thomi $

try
    fclose(self.tmp_file);
    delete('text_test_result.tmp');
catch
end;
