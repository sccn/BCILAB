function [list,indexable] = gui_print_grouplist(grouplist,namemap)
% create a pretty-printed list (cellstr) of grouped data
if ~exist('namemap','var')
    namemap = @(x)x; end

list = {};
% for each group of approaches...
for k=1:2:length(grouplist)
    % put the group in brackets
    list{end+1} = ['[' grouplist{k} ']'];
    sublist = grouplist{k+1};
	for l=1:length(sublist)
        % and indent the sub-items
    	list{end+1} = ['  ' namemap(sublist{l})]; 
    end
end

% create an indexable grouplist, as well...
indexable = grouplist;
indexable(1:2:end) = cellfun(@(x){x},indexable(1:2:end),'UniformOutput',false);
indexable = [indexable{:}];

if nargout == 0
    fprintf('%s\n',list{:}); end