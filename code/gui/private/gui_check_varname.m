function gui_check_varname(edithandle,resetstr)
% check whether the the given edit field contains a valid variable name; warn and reset if not
if ~exist('resetstr','var')
    resetstr = ''; end
if ~isvarname(get(edithandle,'String'))
    errordlg2(sprintf('This field must be a valid variable name.\nOnly letters, numbers and underscores are permitted, and it must begin with a letter.'));
    set(edithandle,'String',resetstr);
end