function second = arg_diff(first,second,use_alternatives,prune_selection)
% Calculate the difference of two arg_specifier struct arrays.
% Difference = arg_diff(First,Second)
%
% In:
%   First : The reference struct array (each element corresponds to one argument, possibly with 
%           non-empty sub-arguments in .children in case of arg_sub*, and possibly with sub-arguments 
%           for alternative choices in .alternatives in case of arg_subswitch and arg_subtoggle)
%
%   Second : A second struct array that is pruned where it is equal to the first (with correct
%            treatment of arg_selection).
%
%   UseAlternatives : whether to consider the contents of the .alternatives field in First when
%                     the value of arg_selection differs in Second (default: true)
%
% Out: 
%   Difference : The second struct array with all entries removed (recursively) that are unchanged
%                from the first.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-10-19

% Copyright (C) Christian Kothe, SCCN, 2013, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

if isempty(first) || isempty(second)
    return; end
if nargin < 3
    use_alternatives = true; end
if nargin < 4
    prune_selection = true; end

% if fields are unequal we need to do a reordered comparison
if length(first) ~= length(second) || ~isequal({first.first_name},{second.first_name})
    remove_first = 1:length(first);
    replicate_second = 1:length(second);
    [dummy,first_in_second] = fast_setdiff({first.first_name},{second.first_name}); %#ok<ASGLU>
    [dummy,second_in_first] = fast_setdiff({second.first_name},{first.first_name}); %#ok<ASGLU>
    remove_first(first_in_second) = [];
    replicate_second(second_in_first) = [];    
    % drop all fields from first that are not in second
    if any(remove_first)
        first(remove_first) = []; end;
    % append the extra fields of second to first, with blank values
    if any(replicate_second)
        dummy = second([]);
        [dummy(1:length(replicate_second)).first_name] = second(replicate_second).first_name;
        first = [first dummy]; 
    end
    % reorder first to match the order of second
    [dummy,idx1] = sort({first.first_name}); %#ok<ASGLU>
    [dummy,idx2] = sort({second.first_name}); %#ok<ASGLU>
    if ~all(idx1 == idx2)
        revidx2(idx2) = 1:length(idx2);
        first = first(idx1(revidx2));
    end
end

% check where the value differs
equal_value = cellfun(@isequalwithequalnans,{first.value},{second.value});
nontrivial = equal_value & ~cellfun('isempty',{first.children});

% at all positions where the value is equal and first has nonempty children, prune the second's children arrays
if any(nontrivial)
    if use_alternatives 
        [second(nontrivial).children] = celldeal(cellfun(@arg_diff,{first(nontrivial).children},{second(nontrivial).children},'UniformOutput',false));
    else
        [second(nontrivial).children] = celldeal(cellfun(@(a,b)arg_diff(a,b,false),{first(nontrivial).children},{second(nontrivial).children},'UniformOutput',false));
    end
end

% at all positions where the value differs, use the corresponding alternative value set (if any)
% as reference
if use_alternatives
    for k=find(~equal_value)
        if ischar(first(k).value) && ~isempty(first(k).alternatives) && length(first(k).alternatives) == length(first(k).range) && any(strcmp(second(k).value,first(k).range))
            second(k).children = arg_diff(first(k).alternatives{strcmp(second(k).value,first(k).range)},second(k).children,use_alternatives,false);
        elseif isequal(first(k).value,false) && length(first(k).alternatives)==2
            second(k).children = arg_diff(first(k).alternatives{1+second(k).value},second(k).children,use_alternatives,false);
        end
    end
end

% determine where the children are empty
empty_children = cellfun('isempty',{second.children});

% if everything is equal, reduce second to its struct fields
% (except if prune_selection is false, in which case we always retain the selection field)
if prune_selection && all(equal_value&empty_children)
    second = second([]);
else
    % otherwise prune everything that's equal except for the selection field
    is_selection = strcmp({second.first_name},'arg_selection');
    second(equal_value & empty_children & ~is_selection) = [];
end
