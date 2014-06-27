function gui_check_scalar(edithandle,resetstr)
% check whether the the given edit field contains a scalar; warn and reset if not
if ~exist('resetstr','var')
    resetstr = ''; end

str = get(edithandle,'String'); num = str2num(str);
if ~isempty(str) && numel(num) ~= 1
    errordlg2('This field must be a number.'); 
    set(edithandle,'String',resetstr);
end