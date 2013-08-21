function reset_texts(self)
%gui_test_result/reset_texts resets the text area.
%
%  Example
%  =======
%  reset_texts is called for example from the constructor of
%  gui_test_result:
%         reset_texts(self);
%
%  See also GUI_TEST_RESULT/GUI_TEST_RESULT.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: reset_texts.m 28 2006-06-11 15:11:59Z thomi $

set(self.text_runs, 'String', ['Runs: 0', ...
    ' / Errors: 0', ...
    ' / Failures: 0']);