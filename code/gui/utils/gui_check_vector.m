function gui_check_vector(edithandle,resetstr)
% check whether the the given edit field contains a row vector; warn and reset if not
if ~exist('resetstr','var')
    resetstr = ''; end

str = get(edithandle,'String');
if ~isempty(str) && ~isvector(str2num(str))
    errordlg2('This field must be of the form [number, number, number, ...]'); 
    set(edithandle,'String',resetstr);
end
