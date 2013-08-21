function self = gui_test_result(progress_bar, text_runs, error_listbox, max_number_of_tests)
%gui_test_result contructor.
%  The constructor creates an object of the class gui_test_result with the
%  following parameters: 
%    progress_bar: handle of the progess bar.
%    text_runs: handle of the text area, which shows the number of test ran,
%               errors and failures.
%    error_listbox: handle of the listbox, which contains the list of 
%                   all errors and failures. 
%    max_number_of_tests: the number of tests.
%
%  Class Info
%  ==========
%  The class gui_test_result is inherited from test_result and shows
%  the test results and the progress at the mlUnit gui. Normally an
%  instance of the class is only created by gui_test_runner, but not by the
%  user.
%
%  Example
%  =======
%  A gui_test_result is created in gui_test_runner/gui:
%         result = gui_test_result(handles.gui_progress_bar, ...
%             handles.gui_text_runs, ...
%             handles.gui_error_list, ...
%             0);
%
%  See also TEST_RESULT, GUI_TEST_RESULT, GUI_TEST_RUNNER/GUI.

self.progress_bar = progress_bar;
self.text_runs = text_runs;
self.error_listbox = error_listbox;
self.max_number_of_tests = max_number_of_tests;
result = test_result;
self = class(self, 'gui_test_result', result);
reset_progress_bar(self);
reset_texts(self);