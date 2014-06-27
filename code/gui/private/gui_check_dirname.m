function gui_check_dirname(edithandle,resetstr)
% check whether the the given edit field contains a valid directory name; warn and reset if not
if ~exist('resetstr','var')
    resetstr = ''; end
if ~exist(get(edithandle,'String'),'dir')
    errordlg2(sprintf('This field must be a valid directory name.'));
    set(edithandle,'String',resetstr);
end