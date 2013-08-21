function reset_progress_bar(self)
%gui_test_result/reset_progress_bar resets the progress bar.
%
%  Example
%  =======
%  reset_progress_bar is called for example from the constructor of
%  gui_test_result:
%         reset_progress_bar(self);
%
%  See also GUI_TEST_RESULT/GUI_TEST_RESULT.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: reset_progress_bar.m 28 2006-06-11 15:11:59Z thomi $

barh(1, 1, 'FaceColor', [1 1 1]);
set(self.progress_bar, 'XLim', [0 1]);
set(self.progress_bar, 'YLim', [0.6 1.4]);
set(self.progress_bar, 'XTick', [], 'XTickLabel', []);
set(self.progress_bar, 'YTick', [], 'YTickLabel', []);
drawnow;
