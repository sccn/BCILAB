function v = testmvnpdf(X,mu,CP)
% Value = logmvnpdf(X,mu,Sigma)
% Derivatives of the log multivariate normal pdf at points X for mean mu and covariance Sigma.
%
% In:
%   X : points at which to evaluate the log-pdf [#dims x #obs]
%
%   mu : mean vector [#dims x 1]
%
%   CP : cholesky factor of precision matrix [#dims x #dims]
%
% Out:
%   Value : the log-pdf evaluated at X [1 x #obs]

%mu = repmat(mu,1,size(X,2)); --> gives an error
F = CP*(X - mu);
v = -0.5*size(X,1)*log(2*pi) + sum(log(diag(CP))) - 0.5*sum(F.*F);

% sample adigator lines:
% d=20; adigator('testmvnpdf',{adigatorCreateDerivInput([d,1],'x'),adigatorCreateDerivInput([d,1],'mu'),adigatorCreateDerivInput([d,d],'CP')},'d_testmvnpdf')
% adigator('d_testmvnpdf',{struct('f',adigatorCreateDerivInput([d,1],'x'),'dx',ones(d,1)),struct('f',adigatorCreateDerivInput([d,1],'mu'),'dmu',ones(d,1)),struct('f',adigatorCreateDerivInput([d,d],'CP'),'dCP',ones(d*d,1))},'dd_testmvnpdf')
%
% res = d_testmvnpdf(struct('f',randn(d,1),'dx',ones(d,1)),struct('f',randn(d,1),'dmu',ones(d,1)),struct('f',randn(d,d),'dCP',ones(d*d,1)))
% res = dd_testmvnpdf(struct('f',randn(d,1),'dx',ones(d,1)),struct('f',randn(d,1),'dmu',ones(d,1)),struct('f',randn(d,d),'dCP',ones(d*d,1)))
