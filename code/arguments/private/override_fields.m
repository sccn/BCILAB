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
    [A_pos,B_pos,B_append] = hlp_nanocache('override',30,@matching,{A.first_name},{B.first_name});
    % save A's children and alternatives
    A_chld = {A(A_pos).children};
    A_alts = {A(A_pos).alternatives};
    % replace A by B's elements
    A(A_pos) = B(B_pos);
    % override the children recursively
    [A(A_pos).children] = celldeal(cellfun(@override_fields,A_chld,{B(B_pos).children},'UniformOutput',false));    
    % override the alternatives recursively
    A_numalts = cellfun('length',A_alts);
    if any(A_numalts)
        B_alts = {B(B_pos).alternatives};
        B_numalts = cellfun('length',B_alts);
        % for each element e in A where there was an alternative
        for e=find(A_numalts)
            % pad the alternatives with [] to ensure same length
            if A_numalts(e) < B_numalts(e)
                [A_alts{(end+1):B_numalts(e)}] = deal([]);
            elseif A_numalts(e) > B_numalts(e)
                [B_alts{(end+1):A_numalts(e)}] = deal([]);
            end
            % override fields at each alternative index
            for a=1:length(A_alts{e})
                A(A_pos(e)).alternatives{a} = override_fields(A_alts{e}{a},B_alts{e}{a}); end
            % [A(A_pos(e)).alternatives] = celldeal(cellfun(@override_fields,A_alts{e},B_alts{e},'UniformOutput',false));
        end
    end    
    % append extra elements of B
    A(end+(1:length(B_append))) = B(B_append);
end


function [A_pos,B_pos,B_append] = matching(A_names,B_names)
% find intersection
[dummy,A_pos,B_pos] = intersect(A_names,B_names); %#ok<ASGLU>
% find part to append
[dummy,B_append] = setdiff(B_names,A_names); %#ok<ASGLU>
