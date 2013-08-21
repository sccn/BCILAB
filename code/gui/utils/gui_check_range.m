function gui_check_range(edithandle,resetstr)
% check whether the the given edit field contains a range; warn and reset if not
if ~exist('resetstr','var')
    resetstr = ''; end

str = get(edithandle,'String');
if ~isempty(str) && numel(str2num(str)) ~= 2
    errordlg2('This field must be of the form [number, number]'); 
    set(edithandle,'String',resetstr);
end
