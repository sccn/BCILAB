% [] = XXXX()
%
% -----------------------------Definition---------------------------------%
%
% description
%
% usage : 
%
% -----------------------------Input--------------------------------------%
%
% XX : description + dimension
%
% -----------------------------Output-------------------------------------%
%
% XX : description + dimension
%
% -----------------------------References---------------------------------%
%
% [1] : XXX
%
%
%   Project : BCI-EEG
%
%   author : A. Barachant
%   date : 2011-XXXX
%   version : 1.0 
%   status : a terminer, termin�   
%   CEA/GRENOBLE-LETI/DTBS
%
%   See also COV, FPCOV, SHCOV, LTCOV, MCDCOV, NORMALIZEDSCM.

% [EOF: XXX.m]

function COV = covariances(X,method_cov,arg_cov)

if (nargin<2)||(isempty(method_cov))
    method_cov = 'scm';
end
if (nargin<3)
    arg_cov = {};
end

[Ne , ~, Nt] = size(X);

COV = zeros(Ne,Ne,Nt);

switch method_cov
    case 'nscm'
        for i=1:Nt
            COV(:,:,i) = NormalizedSCM(X(:,:,i)');
        end
    case 'ntrace'
        for i=1:Nt
            COV(:,:,i) = X(:,:,i)*X(:,:,i)'/trace(X(:,:,i)*X(:,:,i)');
        end
    case 'fp'
        for i=1:Nt
            COV(:,:,i) = FPCov(X(:,:,i)');
        end
    case 'shcov'
        for i=1:Nt
            COV(:,:,i) = shcov(X(:,:,i)');
        end
    case 'shcovft'
        for i=1:Nt
            COV(:,:,i) = shcovft(X(:,:,i)');
        end    
    case 'corr'
        for i=1:Nt
            COV(:,:,i) = corr(X(:,:,i)');
        end
    case 'ltcov'
        for i=1:Nt
            COV(:,:,i) = LTCov2(X(:,:,i),arg_cov{1},arg_cov{2});
            %COV(:,:,i) = LTCov2(X(:,:,i),5,1.5*mean(diag(X(:,:,i)'*X(:,:,i))));
        end
    case 'mcd'
        for i=1:Nt
            res = mcdcov(X(:,:,i)','plots',0);
            COV(:,:,i) = res.cov;
        end
    case 'hann'
        window = repmat(hann(size(X,2)),1,size(X,1));
        for i=1:Nt
            COV(:,:,i) = cov(window.*X(:,:,i)');
        end
    otherwise
        for i=1:Nt
            COV(:,:,i) = cov(X(:,:,i)');
        end
end
