function idx = gui_select_bracketed(listhandle)
% make sure that the user does not select [...] entries (except if not otherwise possible)

contents = cellstr(get(listhandle,'String'));
idx = get(listhandle,'Value');
for k=idx-1:idx-1+length(contents)
    idx = mod(k,length(contents))+1;
    set(listhandle,'Value',idx)
    sel = contents{idx};
    if ~isempty(sel) && sel(1) ~= '['
        break; end
end