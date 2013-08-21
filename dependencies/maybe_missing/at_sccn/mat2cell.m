function c = mat2cell(x,varargin)
%MAT2CELL Break matrix up into a cell array of matrices.
%   C = MAT2CELL(X,M,N) breaks up the 2-D array X into a cell array of  
%   adjacent submatrices of X. X is an array of size [ROW COL], M is the 
%   vector of row sizes (must sum to ROW) and N is the vector of column 
%   sizes (must sum to COL). The elements of M and N determine the size of
%   each cell in C by satisfying the following formula for I = 1:LENGTH(M)
%   and J = 1:LENGTH(N),
%
%       SIZE(C{I,J}) == [M(I) N(J)]
%
%   C = MAT2CELL(X,D1,D2,D3,...,DN) breaks up the multidimensional array X
%   and returns a multidimensional cell array of adjacent submatrices of X.
%   Each of the vector arguments, D1 through DN, should sum to the
%   respective dimension sizes of X, such that for P = 1:N,
%
%       SIZE(X,P) == SUM(DP) 
%
%   The elements of D1 through DN determine the size of each cell in C by
%   satisfying the formula for IP = 1:LENGTH(DP),
%
%       SIZE(C{I1,I2,I3,...,IN}) == [D1(I1) D2(I2) D3(I3) ... DN(IN)]
%
%   C = MAT2CELL(X,R) breaks up an array X by returning a single column
%   cell array, containing the rows of X. R must sum to the number of rows
%   of X. The elements of R determine the size of each cell in C, subject
%   to the following formula for I = 1:LENGTH(R),
%
%       SIZE(C{I},1) == R(I)
%
%   C = MAT2CELL(X,...,[],...) will return an empty cell array whose empty
%   size matches the lengths of the vector arguments, D1 through DN. Note
%   that the length of an empty vector is zero.
%
%   MAT2CELL supports all array types.
%
%	Example:
%	   X = [1 2 3 4; 5 6 7 8; 9 10 11 12];
%	   C = mat2cell(X,[1 2],[1 3])
%	    
%	See also CELL2MAT, NUM2CELL

% Copyright 1984-2007 The MathWorks, Inc.
% $Revision: 1.10.4.5 $  $Date: 2007/10/10 20:51:51 $

if nargin==0
    error('MATLAB:mat2cell:NoInputs',['No input arguments specified. ' ...
            'There should be at least two input arguments.'])
end

switch nargin
    
case 1
    warning('MATLAB:mat2cell:ObsoleteSingleInput', ...
            ['Single input behavior is obsolete and will be removed in a ' ...
            '\n         future release of MATLAB. Use C={X} instead.'])
    c={x};
    
otherwise
    % Allow use of a single input vector argument by assuming the remaining ones
    %   Remaining vectors will be the size of the input array along the
    %   respective dimension.  Necessary for backwards compatibility.
    if nargin == 2
        dims = ndims(x);
        varargin{dims} = 0;
        for n=2:dims
            varargin{n} = size(x,n);
        end
    end
    
    cellsize = cellfun('length',varargin);
    
    % Verify that the number of input arguments meets syntax requirements.
    if length(varargin) ~= ndims(x)
        % Check if the last input vectors have value 1 that would not
        %   affect the output
        vLengthInit = length(varargin);
        singleFlag = 0;
        while (isequal(varargin{end},1)) && (length(varargin) > ndims(x))
            varargin = {varargin{1:end-1}};
            singleFlag = singleFlag+1;
        end
        
        % If trailing ones were found and removed, report this and continue
        if singleFlag
            warning('MATLAB:mat2cell:TrailingUnityVectorArgRemoved', ...
                ['Number of input vectors, %i, did not match the input ' ...
                 'matrix''s number\n         of dimensions, %i. %i ' ...
                 'trailing singleton input vectors were removed.'], ...
                 vLengthInit,ndims(x),singleFlag);
        end
        
        cellsize = cellfun('length',varargin);
        
        if length(varargin) ~= ndims(x)
            error('MATLAB:mat2cell:ArgumentCountMismatch', ...
                ['Number of input vector arguments, %i, does not ' ... 
                    'match the input matrix''s number of dimensions, %i.'], ...
                length(varargin),ndims(x))
        end
    end
    
    % Verify that all dimension size arguments are vectors
    if ~isequal(cellsize,cellfun('prodofsize',varargin))
        error('MATLAB:mat2cell:NonVectorArgument', ...
            'Input arguments, D1 through D%i, should be vectors.', ...
            length(varargin))
    end
    
    % Verify that the input vectors sum to the input matrix dimensions.
    vecsums = zeros(size(varargin));
    for n = 1:length(vecsums)
        vecsums(n) = sum(varargin{n}(:));
    end
    if ~isequal(vecsums,size(x))
        error('MATLAB:mat2cell:VectorSumMismatch', ...
            ['Input arguments, D1 through D%i, must sum to each dimension ' ...
            'of the input matrix size, [%s].'],length(varargin),num2str(size(x)))
    end

    % If an input vector was found to be empty, return an empty cell array
    %   as described in the help. This is necessary for backward
    %   compatibility.
    if any(cellfun('isempty',varargin))
        c = cell(cellfun('length',varargin));
        return
    end

    
    % If matrix is 2-D, execute 2-D code for speed efficiency
    if ndims(x)==2
        rowSizes = varargin{1};
        colSizes = varargin{2};
        rows = length(rowSizes);
        cols = length(colSizes);
        c = cell(rows,cols);
        % Construct each cell element by indexing into X with iterations of 
        %   matrix subscripts (i,j)
        rowStart = 0;
        for i=1:rows
            colStart = 0;
            for j=1:cols
                c{i,j} = x(rowStart+(1:rowSizes(i)),colStart+(1:colSizes(j)));
                colStart = colStart + colSizes(j);
            end
            rowStart = rowStart + rowSizes(i);
        end
        return
    end
    
    % Now treat 3+ dimension arrays
    % Initialize cell array
    c = cell(cellsize);
    
    % Initialize the dimension counter, which keeps track of which cell element 
    %   we are constructing
    dimcounts = ones(size(varargin));
    dimcounts(1) = 0;
    
    % Initialize the subscript references when indexing into the matrix
    %   REFSTART is the set of matrix subscripts to be included in the first 
    %      cell element.
    %   REF is the set of matrix subscripts to be included in the loop's current
    %      cell element.
    %   TREF is a copy of REF, but replaces the 0's with []'s so that the
    %      indexing makes apprpriate empty cells instead of trying to index
    %      with 0, which would cause an error.
    %   REFSTATIC keeps track of which TREF's are [] as the dimension
    %      counter, DIMCOUNTS, increments through each of the cells. This
    %      is required since adding 0 when incrementing the references
    %      returns a value that REF needs to know, but the index should be
    %      [].
    refstart = cell(size(cellsize));
    for n=1:length(refstart)
        refstart{n} = 1:varargin{n}(1);
        if isempty(refstart{n})
            refstart{n} = 0;
        end
    end
    ref = refstart;
    ref{1} = 0;
    tref = ref;
    refstatic = zeros(1,length(cellsize));

    % Construct cell elements by looping through absolute indices and extracting
    %   the appropriate matrix subscripts
    for cind = 1:numel(c)
        % Find the next cell dimension that needs to be incremented
        inc = find(dimcounts<cellsize, 1, 'first');
        % Update the dimension counter using the incremental cell dimension
        dimcounts(1:inc) = [ones(1,inc-1) dimcounts(inc)+1];
 
        % Update the set of matrix subscripts, REF, by indexing into the input 
        %   arguments

        % If not incrementing the first cell dimension, then update REF with the
        %   appropriate earlier cell dimensions as well
        if inc~=1
            [ref{1:inc-1}] = deal(refstart{1:inc-1});
            refstatic(1:inc-1) = 0;
        end
        
        % If we are adding an empty matrix when updating this increment's
        %   reference, then add 0 instead, but index with [], and set the
        %   REFSTATIC tracker.
        if ~varargin{inc}(dimcounts(inc))
            ref{inc} = ref{inc}(end);
            refstatic(inc) = 1;
        else
            % When updating an increment's reference without adding 0,
            %   reset the REFSTATIC tracker
            ref{inc} = ref{inc}(end)+(1:varargin{inc}(dimcounts(inc)));
            refstatic(inc) = 0;
        end
                                
        % If any references are 0, change them to empty so that we are
        % indexing correctly, else let TREF{n} take on REF{n}
        for n = 1:length(refstatic)
            if ~any(ref{n}) || refstatic(n)
                tref{n} = [];
            else
                tref{n} = ref{n};
            end
        end
        
        % Finally construct the cell element by indexing into the matrix
        c{cind} = x(tref{:});

    end
end
