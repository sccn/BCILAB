function A = override_fields(A,B)
% helper function: override specifiers in A selectively by specifiers in B
% Result = override_flags(A,B)
%
% In:
%   A : array of arg_specifier structs
%
%   B : array of arg_specifier structs that should selectively override A
%
% Out:
%   Result : the updated specification
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2014-02-03

% Copyright (C) Christian Kothe, SCCN, 2014, christian@sccn.ucsd.edu
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

if isempty(A)
    A = B;
elseif ~isempty(B)
    % find the positions where B overrides A, and positions in B to append to A
    A_names = {A.first_name};
    B_names = {B.first_name};
    if isequal(A_names,B_names)
        A_pos = 1:length(A_names);
        B_pos = 1:length(B_names);
        B_append = [];
    else
        [A_pos,B_pos,B_append] = hlp_nanocache('override',30,@matching,A_names,B_names);
    end
    % determine whether any element has sub-structure (based on the head)
    A_heads = {A(A_pos).head};
    has_sub = strncmp(A_heads,'arg_sub',7);
    if any(has_sub)
        % initialize A_source with the respective children (ensures right length, etc.)
        [A_source{has_sub}] = A(A_pos(has_sub)).children;
        % for all switches/toggles, we replace the children with the respective alternative
        B_values = {B(B_pos).value};
        is_subtoggle = strcmp(A_heads,'arg_subtoggle');
        is_subswitch = strcmp(A_heads,'arg_subswitch');
        is_selectable = is_subtoggle | is_subswitch;
        if any(is_selectable)
            % get all alternatives
            A_alts = {A(A_pos).alternatives};
            % for each subtoggle that's enabled, use the 'on' alternative as source
            if any(is_subtoggle)
                is_subtoggle = find(is_subtoggle);
                toggle_on = [B_values{is_subtoggle}] ~= 0;
                if any(toggle_on)
                    tmp = vertcat(A_alts{is_subtoggle(toggle_on)});
                    [A_source{is_subtoggle(toggle_on)}] = tmp{:,2};
                end
            end
            % for each subswitch, use the corresponding alternative based on the range
            if any(is_subswitch)
                is_subswitch = find(is_subswitch);
                match = cellfun(@strcmp,B_values(is_subswitch),{A(A_pos(is_subswitch)).range},'UniformOutput',false);
                for m=1:length(match)
                    A_source{is_subswitch(m)} = A_alts{is_subswitch(m)}{match{m}}; end                
            end
        end
        % replace A by B's elements
        A(A_pos) = B(B_pos);
        % merge A_source and B's children recursively
        for k=find(has_sub)
            A(A_pos(k)).children = override_fields(A_source{k},B(B_pos(k)).children); end
        % update alternatives
        if any(is_selectable)
            % write merge results back into the selected alternatives
            if any(is_subtoggle) && any(is_subtoggle(toggle_on))
                for k=find(toggle_on)
                    A_alts{is_subtoggle(k)} = {[],A(A_pos(is_subtoggle(k))).children}; end
            end 
            if any(is_subswitch)
                for m=1:length(match)
                    A_alts{is_subswitch(m)}{match{m}} = A(A_pos(is_subswitch(m))).children; end
            end
            % write back merged alternatives into A
            [A(A_pos(is_selectable)).alternatives] = A_alts{is_selectable};
        end
    else
        % replace A by B's elements
        A(A_pos) = B(B_pos);
    end
    % append extra elements of B
    A(end+(1:length(B_append))) = B(B_append);
end


function [A_pos,B_pos,B_append] = matching(A_names,B_names)
% find intersection
[dummy,A_pos,B_pos] = intersect(A_names,B_names); %#ok<ASGLU>
% find part to append
[dummy,B_append] = setdiff(B_names,A_names); %#ok<ASGLU>
