function x = kth_sequence(x,k)
% Generate a rearranged sequence of items to process for the k'th participating process (counting from 0).
%
% Each sequence covers all items, but they are chosen to first cover items that are maximally far away
% from those of previous sequences (e.g., {forward,reversed,inside-out,from-center-left,from-center-right, ...})

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

x = x(center_out_sequence(length(x),nested_index(length(x),k)));

function pos = nested_index(len,k)
% Generate the k'th index of a nested sequence of indices in the 
% sequence 1:len (1,len,len/2,1/4*len,3/4*len, ...).
% This function is not perfect (i.e., there can be duplicates) but any such issue is 
% handled gracefully in the context where it is used (with overlapping center_out_sequence's).
if k == 0
    pos = 1;
elseif k == 1
    pos = len;
else    
    num_covered = 2;                           % number of indices already covered in previous tiers
    for tier = 1:50                            % these are nesting tiers (like rows in a balanced tree)
        num_indices = 2^(tier-1);              % number of new indices in this tier
        if k<num_covered+num_indices           % if k is in this tier
            numerator = (k-num_covered+1)*2-1; % numerator of the k'th fractional position in this tier
            denominator = 2^tier;              % denominator of the k'th fractional position
            pos = ceil(len*numerator/denominator);
            break;
        end
        num_covered = num_covered+num_indices;
    end
end

pos = min(max(pos,1),len);
    
function seq = center_out_sequence(len,ctr)
% Generate a center-out sequence beginning at ctr that covers the range 1:len.
seq = ctr;
for offset=1:len
    for side=[-1,+1]
        pos = ctr + offset*side;
        if pos>=1 && pos <=len && ~ismember(pos,seq)
            seq(end+1) = pos; end
        if length(seq) == len
            return; end
    end
end
