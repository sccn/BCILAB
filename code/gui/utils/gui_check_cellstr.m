function gui_check_cellstr(edithandle,resetstr)
% check whether the the given edit field contains a valid cell-string array; warn and reset if not
if ~exist('resetstr','var')
    resetstr = ''; end
try
    data = evalin('base',get(edithandle,'String'));
    assert(iscellstr(data),'Error!');
catch
    errordlg2(sprintf('This field must be a cell array of strings.\nExample: {''string1'',''string2'',''string3'', ...}'));
    set(edithandle,'String',resetstr);    
end
