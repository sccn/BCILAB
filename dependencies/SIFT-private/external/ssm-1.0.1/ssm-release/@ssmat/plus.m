function C = plus(A, B)

%@SSMAT/PLUS Addition of SSMAT objects.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if ~isa(A, 'ssmat'), if isnumeric(A) && ndims(A) <= 3, A = ssmat(A); else error('ssm:ssmat:plus:UnableToConvert', 'Input cannot be converted to SSMAT class.'); end
elseif ~isa(B, 'ssmat'), if isnumeric(B) && ndims(B) <= 3, B = ssmat(B); else error('ssm:ssmat:plus:UnableToConvert', 'Input cannot be converted to SSMAT class.'); end
end
if size(A) ~= size(B), error('ssm:ssmat:plus:dimagree', 'State space matrix dimensions must agree.'); end

if isempty(A.mmask), Ammask = false(size(A.mat)); else Ammask = A.mmask; end
if isempty(B.mmask), Bmmask = false(size(B.mat)); else Bmmask = B.mmask; end

if isempty(A.dmmask) && isempty(B.dmmask)
    C.mat       = A.mat + B.mat;
    C.mmask     = Ammask | Bmmask;
    C.dmmask    = [];
    C.dvec      = zeros(0, 1);
    C.dvmask    = [];
elseif ~isempty(A.dmmask) && ~isempty(B.dmmask)
    if size(A.dvec, 2) ~= size(B.dvec, 2), error('ssm:ssmat:plus:dimagree', 'State space matrix dimensions must agree.'); end
    n                               = size(A.dvec, 2);
    Amat                            = repmat(A.mat, [1 1 n]);
    Amat(repmat(A.dmmask, [1 1 n])) = A.dvec;
    Bmat                            = repmat(B.mat, [1 1 n]);
    Bmat(repmat(B.dmmask, [1 1 n])) = B.dvec;
    Cmat                            = Amat + Bmat;
    C.mat                           = Cmat(:,:, 1);
    C.mat(A.dmmask | B.dmmask)      = 0;
    C.mmask                         = Ammask | Bmmask;
    C.dmmask                        = A.dmmask | B.dmmask;
    C.dvec                          = reshape(Cmat(repmat(C.dmmask, [1 1 n])), nnz(C.dmmask), n);
    if isempty(A.dvmask) && isempty(B.dvmask)
        C.dvmask            = [];
    else
        if isempty(A.dvmask), A.dvmask = false(size(A.dvec, 1), 1); end
        if isempty(B.dvmask), B.dvmask = false(size(B.dvec, 1), 1); end
        Advmask             = A.dmmask;
        Advmask(Advmask)    = A.dvmask;
        Bdvmask             = B.dmmask;
        Bdvmask(Bdvmask)    = B.dvmask;
        Cdvmask             = Advmask | Bdvmask;
        C.dvmask            = Cdvmask(C.dmmask);
    end
else
    if isempty(B.dmmask)
        temp    = B;
        B       = A;
        A       = temp;
    end
    n                               = size(B.dvec, 2);
    Amat                            = repmat(A.mat, [1 1 n]);
    Bmat                            = repmat(B.mat, [1 1 n]);
    Bmat(repmat(B.dmmask, [1 1 n])) = B.dvec;
    Cmat                            = Amat + Bmat;
    C.mat                           = Cmat(:,:, 1);
    C.mat(B.dmmask)                 = 0;
    C.mmask                         = Ammask | Bmmask;
    C.dmmask                        = B.dmmask;
    C.dvec                          = reshape(Cmat(repmat(C.dmmask, [1 1 n])), size(B.dvec));
    C.dvmask                        = B.dvmask;
end

%% Register object instance %%
C   = class(C, 'ssmat');

