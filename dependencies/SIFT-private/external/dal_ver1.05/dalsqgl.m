% dalsqgl - DAL with squared loss and grouped L1 regularization
%
% Overview:
%  Solves the optimization problem:
%   xx = argmin 0.5||A*x-bb||^2 + lambda*||x||_G1
%  where
%   ||x||_G1 = sum(sqrt(sum(xx.^2)))
%  (grouped L1 norm)
%
% Syntax:
%  [xx,status]=dalsqgl(xx0, A, bb, lambda, <opt>)
%
% Inputs:
%  xx0    : initial solution ([nn,1] with opt.blks or [ns nc] with
%           ns*nc=nn for nc groups of size ns)
%  A      : the design matrix A ([mm,nn] or [mm,ns,nc])
%  bb     : the target vector ([mm,1])
%  lambda : the regularization constant
%  <opt>  : list of 'fieldname1', value1, 'filedname2', value2, ...
%   blks     : vector that contains the size of the groups. 
%              sum(opt.blks)=nn. If omitted, opt.blks = [ns,..., ns]
%              and length(opt.blks)=nc, where nc is the number of groups.
%   stopcond : stopping condition, which can be
%              'pdg'  : Use relative primal dual gap (default)
%              'fval' : Use the objective function value
%           (see dal.m for other options)
% Outputs:
%  xx     : the final solution ([ns,nc] or [nn,1])
%  status : various status values
%
% Example:
% m = 1024; n = [64 64]; k = round(0.2*n(1)); A=randn(m,prod(n));
% w0=randsparse(n,k); bb=A*w0(:)+0.01*randn(m,1);
% lambda=0.1*max(sqrt(sum(reshape(A'*bb,n).^2)));
% [ww,stat]=dalsqgl(zeros(n), A, bb, lambda);
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [ww,status]=dalsqgl(ww, A, bb, lambda, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt,'solver','cg',...
                     'stopcond','pdg',...
                     'blks',[]);

if isempty(opt.blks)
  opt.blks=size(ww,1)*ones(1,size(ww,2));
  ww = ww(:);
end

prob.floss    = struct('p',@loss_sqp,'d',@loss_sqd,'args',{{bb}});
prob.fspec    = @gl_spec;
prob.dnorm    = @gl_dnorm;
prob.obj      = @objdalgl;
prob.softth   = @gl_softth;
prob.stopcond = opt.stopcond;
prob.ll       = -inf*ones(size(bb));
prob.uu       = inf*ones(size(bb));
prob.Ac       =[];
prob.bc       =[];
prob.info     = struct('blks',opt.blks);

if isequal(opt.solver,'cg')
  prob.hessMult = @hessMultdalgl;
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

if all(opt.blks==opt.blks(1))
  ns=opt.blks(1);
  nc=length(ww)/ns;
  ww=reshape(ww, [ns,nc]);
end
