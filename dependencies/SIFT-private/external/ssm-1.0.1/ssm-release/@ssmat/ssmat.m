function matrix = ssmat(m, mmask, dmmask, dvec, dvmask)

%@SSMAT/SSMAT State space matrix class constructor.
%   matrix = SSMAT(m[, mmask]) creates a stationary or dynamic matrix.
%   matrix = SSMAT(m, mmask, dmmask[, dvec, dvmask]) creates a dynamic matrix.
%       m is a 2-d or 3-d matrix, where the third dimension is time. If 2-d
%           then m is the constant and stationary part of matrix.
%       mmask is a 2-d logical matrix masking the variable part of matrix.
%           (Set to [] for constant state space matrices)
%       dmmask is a logical matrix masking the dynamic part of matrix.
%       dvec is a nnz(dmmask)*n matrix, dvec(:, t) is the values for the
%           dynamic part at time t, default is nnz(dmmask)*1 zero matrix.
%           (m(dmmask) = dvec(:, t) for time t)
%       dvmask is a nnz(dmmask)*1 logical vector masking the variable part of
%           the dynamic vector.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin == 0
    %% Default constructor %%
    matrix.mat      = [];
    matrix.mmask    = [];
    matrix.dmmask   = [];
    matrix.dvec     = zeros(0, 1);
    matrix.dvmask   = [];
elseif isa(m, 'ssmat')
    %% Copy constructor %%
    matrix          = m;
    return;
else
    if ~isnumeric(m) || ndims(m) > 3, error('ssm:ssmat:ssmat:InputError', 'm must be a 2-d or 3-d matrix.'); end
    if nargin < 2 || isempty(mmask), mmask = []; elseif ~islogical(mmask) || ~isequal([size(m, 1) size(m, 2)], size(mmask)), error('ssm:ssmat:ssmat:InputError', 'mmask must be a logical matrix with the same size as m.'); end

    if ndims(m) == 3
        %% Dynamic SSMAT (will ignore the other arguments) %%
        matrix.mat      = zeros(size(m, 1), size(m, 2));
        matrix.mmask    = mmask;
        matrix.dmmask   = true(size(m, 1), size(m, 2));
        matrix.dvec     = reshape(m, size(m, 1)*size(m, 2), size(m, 3));
        matrix.dvmask   = [];
    elseif nargin < 3
        %% Stationary SSMAT %%
        matrix.mat      = m;
        matrix.mmask    = mmask;
        matrix.dmmask   = [];
        matrix.dvec     = zeros(0, 1);
        matrix.dvmask   = [];
    else
        if isempty(dmmask), dmmask = []; elseif ~islogical(dmmask) || ~isequal(size(m), size(dmmask)), error('ssm:ssmat:ssmat:InputError', 'dmmask must be a logical matrix with the same size as m.'); end
        if ~isempty(mmask) && ~isempty(dmmask) && any(any(mmask & dmmask)), error('ssm:ssmat:ssmat:InputError', 'mmask and dmmask cannot overlap (use dvmask for variable dynamic elements).'); end
        if nargin < 4 || isempty(dvec), dvec = zeros(nnz(dmmask), 1); elseif ~isnumeric(dvec) || size(dvec, 1) ~= nnz(dmmask), error('ssm:ssmat:ssmat:InputError', 'dvec must be a non-empty matrix with nnz(dmmask) rows.'); end
        if nargin < 5 || isempty(dvmask), dvmask = []; elseif ~islogical(dvmask) || ~isequal(size(dvmask), [size(dvec, 1) 1]), error('ssm:ssmat:ssmat:InputError', 'dvmask must be a logical vector of size nnz(dmmask)*1.'); end
        
        %% Dynamic SSMAT %%
        matrix.mat      = m;
        matrix.mmask    = mmask;
        matrix.dmmask   = dmmask;
        matrix.dvec     = dvec;
        matrix.dvmask   = dvmask;
    end
end

%% Eliminate mask degeneracy %%
if ~isempty(matrix.mmask) && all(~matrix.mmask(:)), matrix.mmask = []; end
if ~isempty(matrix.dmmask) && all(~matrix.dmmask(:)), matrix.dmmask = []; end
if ~isempty(matrix.dvmask) && all(~matrix.dvmask(:)), matrix.dvmask = []; end

%% Register object instance %%
matrix  = class(matrix, 'ssmat');

