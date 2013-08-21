function progress_bar(self)
%gui_test_result/progres_bar updates the progress bar of the mlUnit gui.
%  The length of the progress bar is defined through the number of tests 
%  and with each executed tests the progress is increased.
%
%  Example
%  =======
%  progress_bar is called for example from gui_test_result/update:
%         progress_bar(self);
%
%  See also GUI_TEST_RESULT/UPDATE.

%  This Software and all associated files are released unter the 
%  GNU General Public License (GPL), see LICENSE for details.
%  
%  §Author: Thomas Dohmke <thomas@dohmke.de> §
%  $Id: progress_bar.m 83 2006-10-10 19:06:10Z thomi $

runs = get_tests_run(self);
axes(self.progress_bar);
if ((get_errors(self) > 0) || (get_failures(self) > 0))
    barh(1, runs, 'FaceColor', [1 0 0]);
else
    barh(1, runs, 'FaceColor', [0 1 0]);
end;

set(self.progress_bar, 'XLim', [0 self.max_number_of_tests]);
set(self.progress_bar, 'YLim', [0.6 1.4]);
set(self.progress_bar, 'XTick', [], 'XTickLabel', []);
set(self.progress_bar, 'YTick', [], 'YTickLabel', []);
