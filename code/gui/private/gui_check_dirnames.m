function gui_check_dirnames(edithandle,resetstr)
% check whether the the given edit field contains a valid directory name; warn and reset if not
if ~exist('resetstr','var')
    resetstr = ''; end
try
    str = get(edithandle,'String');
    if exist(str,'dir')
        set(edithandle,'String',hlp_tostring({str}));
    else
        data = evalin('base',str);
        assert(iscellstr(data),'Error!');
        assert(all(cellfun(@(x)exist(x,'dir'),data)));
    end
catch
    errordlg2('This field must be a valid directory name or cell array of directory names.');
    set(edithandle,'String',resetstr);    
end
