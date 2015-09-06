% dalsqen - DAL with squared loss and the Elastic-net regularization
%
% Overview:
%  Solves the optimization problem:
%   [xx, bias] = argmin 0.5||A*x-bb||^2 + lambda*sum(theta*abs(x)+0.5*(1-theta)*x.^2)
%
% Syntax:
%  [xx,status]=dalsqen(xx, A, bb, lambda, theta, <opt>)
%
% Inputs:
%  xx     : initial solution ([nn,1])
%  A      : the design matrix A ([mm,nn]) or a cell array {fA, fAT, mm, nn}
%           where fA and fAT are function handles to the functions that
%           return A*x and A'*x, respectively, and mm and nn are the
%           numbers of rows and columns of A.
%  bb     : the target vector ([mm,1])
%  lambda : the regularization constant
%  theta  : parameter controlling the balance between the
%           L1 term and the L2 term (theta=0->L2, theta=1->L1)
%  <opt>  : list of 'fieldname1', value1, 'filedname2', value2, ...
%   stopcond : stopping condition, which can be
%              'pdg'  : Use relative primal dual gap (default)
%              'fval' : Use the objective function value
%           (see dal.m for other options)
% Outputs:
%  xx     : the final solution ([nn,1])
%  status : various status values
%
% Example:
% m = 1024; n = 4096; k = round(0.04*n); A=randn(m,n);
% w0=randsparse(n,k); bb=A*w0+0.01*randn(m,1);
% lambda=0.1*max(abs(A'*bb));
% [ww,stat]=dalsqen(zeros(n,1), A, bb, lambda, 0.5);
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [ww,status]=dalsqen(ww, A, bb, lambda, theta, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt,'solver','nt',...
                     'stopcond','pdg');



prob.floss    = struct('p',@loss_sqp,'d',@loss_sqd,'args',{{bb}});
prob.fspec    = @(ww)en_spec(ww, theta);
prob.dnorm    = @(vv)en_dnorm(vv, lambda, theta);
prob.obj      = @objdalen;
prob.softth   = @en_softth;
prob.stopcond = opt.stopcond;
prob.ll       = -inf*ones(size(bb));
prob.uu       = inf*ones(size(bb));
prob.Ac       =[];
prob.bc       =[];
prob.info     =struct('theta',theta);

if isequal(opt.solver,'cg')
  prob.hessMult = @hessMultdalen;
end

if isequal(opt.stopcond,'fval')
  opt.feval = 1;
end

if isnumeric(A)
  A = A(:,:);
  [mm,nn]=size(A);
  At=A';
  fA = struct('times',@(x)A*x,...
              'Ttimes',@(x)At*x,...
              'slice', @(I)A(:,I));
  clear At;
elseif iscell(A)
  mm = A{3};
  nn = A{4};
  fAslice = @(I)fA(sparse(I,1:length(I),ones(length(I),1), nn, length(I)));
  fA = struct('times',A{1},...
              'Ttimes',A{2},...
              'slice',fAslice);
else
  error('A must be either numeric or cell {@(x)A*x, @(y)(A''*y), mm, nn}');
end

prob.mm       = mm;
prob.nn       = nn;

[ww,uu,status]=dal(prob,ww,[],fA,[],lambda,opt);

