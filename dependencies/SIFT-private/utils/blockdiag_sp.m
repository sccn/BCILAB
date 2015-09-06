function d = blockdiag_sp( varargin )
%BLOCKDIAG Build a block diagonal matrix of the input arguments.
%
% BLOCKDIAG( A, B, C, ... ) returns a block diagonal matrix which has
% the matrices A, B, C, ... on the main diagonal.
%
% BLOCKDIAG( A, B, C, 'sparse' ) returns the same, but as a sparse
% matrix.
%
% See also DIAG, BLKDIAG, HORZCAT, VERTCAT.

% Author: Peter J. Acklam
% Time-stamp: 1999-03-10 17:17:27
% E-mail: jacklam@math.uio.no
% WWW URL: http://www.math.uio.no/~jacklam

   if length( varargin ) < 1
      error( 'Not enough input arguments.' );
   end

   % Should we return a sparse matrix?
   want_sparse = strcmp( varargin{end}, 'sparse' );
   if want_sparse, varargin = varargin(1:end-1); end

   nargsin = length( varargin );
   if nargsin < 1
      error( 'Not enough input arguments.' );
   end

   % Calculate the size of the output matrix.
   rows = 0; % Number of rows in output matrix.
   cols = 0; % Number of columns in output matrix.
   nnz = 0; % Number of non-zero elements in output matrix.
   for k = 1:nargsin
      [ r, c, d ] = size( varargin{k} );
      if d > 1
         error( 'Matrices can not have higher dimesion than 2' );
      end
      rows = rows + r;
      cols = cols + c;
      if want_sparse
         nnz = nnz + length( find( varargin{k} ) );
         varargin{k} = sparse( varargin{k} );
      end
   end

   % Initialize output matrix.
   if want_sparse
      d = sparse( [], [], [], rows, cols, nnz );
      if nnz == 0, return, end
   else
      d = zeros( rows, cols );
   end

   % Fill the input matrices into the output matrix.
   i = 0;
   j = 0;
   for k = 1:nargsin
      [ r, c ] = size( varargin{k} );
      d( i+1:i+r , j+1:j+c ) = varargin{k};
      i = i + r;
      j = j + c;
   end