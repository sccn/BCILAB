function fp = hlp_fingerprint(data)
% Make a fingerprint (hash) of the given data structure.
% Fingerprint = hlp_fingerprint(Data)
%
% This includes all contents; however, large arrays (such as EEG.data) are only spot-checked. For
% thorough checking, use hlp_cryptohash.
%
% In:
%   Data        : some data structure
%
% Out:
%   Fingerprint : an integer that identifies the data
%
% Notes:
%   The fingerprint is not unique and identifies the data set only with a certain (albeit high)
%   probability. 
%
%   On MATLAB versions prior to 2008b, hlp_fingerprint cannot be used concurrently from timers,
%   and also may alter the random generator's state if cancelled via Ctrl+C.
%
% Examples:
%   % calculate the hash of a large data structure
%   hash = hlp_fingerprint(data);
%
% See also:
%   hlp_cryptohash
%
%                                   Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                   2010-04-02

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
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

warning off MATLAB:structOnObject

if hlp_matlab_version >= 707
    fp = fingerprint(data,RandStream('swb2712','Seed',5183));
else
    try
        % save & override random state
        randstate = rand('state'); %#ok<*RAND>
        rand('state',5183);
        % make the fingerprint
        fp = fingerprint(data,0);
        % restore random state
        rand('state',randstate);
    catch e
        % restore random state in case of an error...
        rand('state',randstate);
        rethrow(e);
    end
end


% make a fingerprint of the given data structure
function fp = fingerprint(data,rs)
% convert data into a string representation
data = summarize(data,rs);
% make sure that it does not contain 0's
data(data==0) = 'x';
% obtain a hash code via Java (MATLAB does not support proper integer arithmetic...)
str = java.lang.String(data);
fp = str.hashCode()+2^31;


% get a recursive string summary of arbitrary data
function x = summarize(x,rs)
if isnumeric(x)
    % numeric array
    if ~isreal(x)
        x = [real(x) imag(x)]; end
    if issparse(x)
        try
            x = [find(x) nonzeros(x)];
        catch
            x = [find(x)' nonzeros(x)];
        end
    end
    if numel(x) <= 4096
        % small matrices are hashed completely
        try
            x = ['n' typecast([size(x) x(:)'],'uint8')];
        catch
            if hlp_matlab_version <= 702
                x = ['n' typecast([size(x) double(x(:))'],'uint8')]; end
        end
    else
        % large matrices are spot-checked
        ne = numel(x);
        count = floor(256 + (ne-256)/1000);
        if hlp_matlab_version < 707
            indices = 1+floor((ne-1)*rand(1,count));
        else
            indices = 1+floor((ne-1)*rand(rs,1,count));
        end
        if size(x,2) == 1
            % x is a column vector: reindexed expression needs to be transposed
            x = ['n' typecast([size(x) x(indices)'],'uint8')];
        else
            % x is a matrix or row vector: shape follows that of indices
            x = ['n' typecast([size(x) x(indices)],'uint8')];
        end
    end
elseif iscell(x)
    % cell array
    sizeprod = cellfun('prodofsize',x(:));
    if all(sizeprod <= 1) && any(sizeprod)
        % all scalar elements (some empty, but not all)
        if all(cellfun('isclass',x(:),'double')) || all(cellfun('isclass',x(:),'single'))
            % standard floating-point scalars
            if cellfun('isreal',x)
                % all real
                x = ['cdr' typecast([size(x) x{:}],'uint8')];
            else
                % some complex
                x = ['cdc' typecast([size(x) real([x{:}]) imag([x{:}])],'uint8')];
            end
        elseif cellfun('isclass',x(:),'logical')
            % all logical
            x = ['cl' typecast(uint32(size(x)),'uint8') uint8([x{:}])];
        elseif cellfun('isclass',x(:),'char')
            % all single chars
            x = ['ccs' typecast(uint32(size(x)),'uint8') x{:}];
        else
            % generic types (structs, cells, integers, handles, ...)
            tmp = cellfun(@summarize,x,repmat({rs},size(x)),'UniformOutput',false);
            x = ['cg' typecast(uint32(size(x)),'uint8') tmp{:}];
        end
    elseif isempty(x)
        % empty cell array
        x = ['ce' typecast(uint32(size(x)),'uint8')];
    else
        % some non-scalar elements
        dims = cellfun('ndims',x(:));
        size1 = cellfun('size',x(:),1);
        size2 = cellfun('size',x(:),2);
        if all((size1+size2 == 0) & (dims == 2))
            % all empty and nondegenerate elements
            if all(cellfun('isclass',x(:),'double'))
                % []'s
                x = ['ced' typecast(uint32(size(x)),'uint8')];
            elseif all(cellfun('isclass',x(:),'cell'))
                % {}'s
                x = ['cec' typecast(uint32(size(x)),'uint8')];
            elseif all(cellfun('isclass',x(:),'struct'))
                % struct()'s
                x = ['ces' typecast(uint32(size(x)),'uint8')];
            elseif length(unique(cellfun(@class,x(:),'UniformOutput',false))) == 1
                % same class
                x = ['cex' class(x{1}) typecast(uint32(size(x)),'uint8')];
            else
                % arbitrary class...
                tmp = cellfun(@summarize,x,repmat({rs},size(x)),'UniformOutput',false);
                x = ['cg' typecast(uint32(size(x)),'uint8') tmp{:}];
            end
        elseif all((cellfun('isclass',x(:),'char') & size1 <= 1) | (sizeprod==0 & cellfun('isclass',x(:),'double')))
            % all horizontal strings or proper empty strings, possibly some []'s
            x = ['cch' [x{:}] typecast(uint32(size2'),'uint8')];
        else
            % arbitrary sizes...
            if all(cellfun('isclass',x(:),'double')) || all(cellfun('isclass',x(:),'single'))
                % all standard floating-point types...
                tmp = cellfun(@vectorize,x,'UniformOutput',false);
                % treat as a big vector...
                x = ['cn' typecast(uint32(size(x)),'uint8')  summarize([tmp{:}],rs)];
            else            
                tmp = cellfun(@summarize,x,repmat({rs},size(x)),'UniformOutput',false);
                x = ['cg' typecast(uint32(size(x)),'uint8') tmp{:}];
            end
        end
    end
elseif ischar(x)
    % char array
    x = ['c' x(:)'];    
elseif isstruct(x)
    % struct
    fn = fieldnames(x)';
    if numel(x) > length(fn)
        % summarize over struct fields to expose homogeneity
        x = cellfun(@(f)summarize({x.(f)},rs),fn,'UniformOutput',false);
        x = ['s' [fn{:}] ':' [x{:}]];
    else
        % summarize over struct elements
        x = ['s' [fn{:}] ':' summarize(struct2cell(x),rs)];
    end
elseif islogical(x)
    % logical array
    x = ['l' typecast(uint32(size(x)),'uint8') uint8(x(:)')];
elseif isa(x,'function_handle')
    if strncmp(char(x),'@(',2)
        f = functions(x);
        x = ['f' f.function summarize(f.workspace{1},rs)];
    else
        x = ['f ' char(x)];
    end
elseif isobject(x)
    x = ['o' class(x) ':' summarize(struct(x),rs)];
else
    try
        x = ['u' class(x) ':' summarize(struct(x),rs)];
    catch
        warning('BCILAB:hlp_fingerprint:unsupported_type','Unsupported type: %s',class(x));
        error; %#ok<LTARG>
    end
end


function x = vectorize(x)
x = x(:)';
