function gui_check_fraction(edithandle,resetstr)
% check whether the the given edit field contains a fraction (or percentage); warn and reset if not
if ~exist('resetstr','var')
    resetstr = ''; end

str = get(edithandle,'String'); 
if ~isempty(str) && str(end) == '%'
    num = str2num(str(1:end-1))/100;
else
    num = str2num(str);
end
if ~isempty(str) && (numel(num) ~= 1 || num < 0 || num > 1)
    errordlg2('This field must be a fraction or percentage.'); 
    set(edithandle,'String',resetstr);
end