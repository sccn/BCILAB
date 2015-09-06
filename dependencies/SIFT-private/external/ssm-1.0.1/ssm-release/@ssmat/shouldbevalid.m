function valid = shouldbevalid(A)

%@SSMAT/SHOULDBEVALID State space matrix should be valid.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

mat_type        = isnumeric(A.mat) && ndims(A.mat) == 2;
mmask_type      = isempty(A.mmask) || (islogical(A.mmask) && isequal(size(A.mmask), size(A.mat)));
dmmask_type     = isempty(A.dmmask) || (islogical(A.dmmask) && isequal(size(A.dmmask), size(A.mat)));
dvec_type       = isnumeric(A.dvec) && ndims(A.dvec) == 2 && size(A.dvec, 1) == nnz(A.dmmask) && size(A.dvec, 2) > 0;
dvmask_type     = isempty(A.dvmask) || (islogical(A.dvmask) && isequal(size(A.dvmask), [nnz(A.dmmask) 1]));

mmask_content   = isempty(A.mmask) || any(A.mmask(:));
dmmask_content  = isempty(A.dmmask) || any(A.dmmask(:));
dvmask_content  = isempty(A.dvmask) || any(A.dvmask);

consistency     = isempty(A.mmask) || isempty(A.dmmask) || ~any(any(A.mmask & A.dmmask));

type            = mat_type && mmask_type && dmmask_type && dvec_type && dvmask_type;
content         = mmask_content && dmmask_content && dvmask_content;
valid           = type && content && consistency;

if ~valid, warning('ssm:ssmat:shouldbevalid', 'Invalid state space matrix.'); end

