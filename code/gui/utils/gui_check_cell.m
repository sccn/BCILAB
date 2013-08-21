function gui_check_cell(edithandle,resetstr)
% check whether the the given edit field contains a valid cell array; warn and reset if not
if ~exist('resetstr','var')
    resetstr = ''; end
try
    data = evalin('base',get(edithandle,'String'));
    assert(iscell(data),'Error!');
catch
    errordlg2(sprintf('This field must be a cell array.\nExample: {''test'',3,[4,5], ...}'));
    set(edithandle,'String',resetstr);    
end
