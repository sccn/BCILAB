function r = isequal_weak(a,b)
% Compare 2 values for equality under relaxed conditions.
% Equal = isequal_weak(A,B)
%
% In:
%   A : first value
%
%   B : second value
%
% Out:
%   Equal : whether the two values are weakly equal, ignoring:
%           * NaNs
%           * contents of Java objects
%           * variables accessible to anonymous functions that are not referenced in the code
%           * new-style classes that have the same struct() representation
%
% Examples:
%   % test two data sets for equivalence
%   isequal_weak(EEG1,EEG2)
%
% See also:
%   isequal, isequalwithequalnans
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-02-19

% Copyright (C) Christian Kothe, SCCN, 2011, christian@sccn.ucsd.edu
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

if isequal(a,b)
    r = true;
else
    ca = class(a);
    if ~strcmp(ca,class(b))
        r = false;
        % unequal class
        return;
    else
        % same class; check in detail
        switch ca
            case {'double','char','logical'}
                % common non-recursive data structure
                r = isequalwithequalnans(a,b);
            case 'cell'
                % recurse into cells
                if ~isequal(size(a),size(b))
                    r = false;
                else
                    for k=1:numel(a)
                        if ~isequal_weak(a{k},b{k})
                            r = false;
                            return;
                        end
                    end
                    r = true;
                end
            case 'struct'
                % recurse into structs
                if ~isequal(size(a),size(b)) || ~isequal(fieldnames(a),fieldnames(b))
                    r = false;
                else
                    r = isequal_weak(struct2cell(a),struct2cell(b));
                end
            case 'function_handle'
                if isequal(a,b)
                    r = true;
                else
                    % functions are considered different by isequal...
                    sa = char(a);
                    if strcmp(sa,char(b)) && sa(1) == '@'
                        % but are both anonymous functions with identical code:
                        fa = functions(a);
                        fb = functions(b);
                        % equality is determined by referenced variables
                        if ~isempty(fa.workspace) && ~isempty(fb.workspace)
                            r = isequal_weak(fa.workspace{1},fb.workspace{1});
                        elseif isempty(fa.workspace) && isempty(fb.workspace)
                            r = true;
                        elseif (isempty(fa.workspace) && isempty(fieldnames(fb.workspace{1}))) || (isempty(fb.workspace) && isempty(fieldnames(fa.workspace{1})))
                            r = true;
                        else
                            r = false;
                        end
                    else
                        r = false;
                    end
                end
            otherwise
                % misc class
                if isnumeric(a)
                    % misc number format
                    r = isequalwithequalnans(a,b);
                elseif isobject(a)
                    % MATLAB objects with same class: compare as structs
                    r = isequal_weak(struct(a),struct(b));
                elseif isjava(a)
                    % Java objects with same class: ignore contents
                    r = true;
                else                
                    try
                        % new-style classes with the same class?
                        r = isequal_weak(struct(a),struct(b));
                    catch
                        % some unknown data type
                        r = isequalwithequalnans(a,b);
                    end
                end
        end
    end
end
