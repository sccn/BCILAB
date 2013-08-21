% dallrgl - DAL with logistic loss and grouped L1 regularization
%
% Overview:
%  Solves the optimization problem:
%   [xx,bias] = argmin sum(log(1+exp(-yy.*(A*x+bias)))) + lambda*||x||_G1
%  where
%   ||x||_G1 = sum(sqrt(sum(xx(Ii).^2)))
%   (Ii is the index-set of the i-th group
%
% Syntax:
%  [xx,bias,status]=dallrgl(xx0, bias0, A, yy, lambda, <opt>)
%
% Inputs:
%  xx0    : initial solution ([nn,1] with opt.blks or [ns nc] with
%           ns*nc=nn for nc groups of size ns)
%  bias0  : initial bias (set [] if bias term is unnecessary)
%  A      : the design matrix A ([mm,nn]) or a cell array {fA, fAT, mm, nn}
%           where fA and fAT are function handles to the functions that
%           return A*x and A'*x, respectively, and mm and nn are the
%           numbers of rows and columns of A.
%  yy     : the target label vector (-1 or +1) ([mm,1])
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
%  xx     : the final solution ([nn,1])
%  bias   : the final bias term (scalar)
%  status : various status values
%
% Example:
% m = 1024; n = [64 64]; k = round(0.1*n(1)); A=randn(m,prod(n));
% w0=randsparse(n,k); yy=sign(A*w0(:)+0.01*randn(m,1));
% lambda=0.1*max(sqrt(sum(reshape(A'*yy/2,n).^2)));
% [ww,bias,stat]=dallrgl(zeros(n), 0, A, yy, lambda);
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [ww,bias,status]=dallrgl(ww, bias, A, yy, lambda, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt,'solver','cg',...
                     'stopcond','pdg',...
                     'blks',[]);

if isempty(opt.blks)
  opt.blks=size(ww,1)*ones(1,size(ww,2));
  ww = ww(:);
end

prob.floss    = struct('p',@loss_lrp,'d',@loss_lrd,'args',{{yy}});
prob.fspec    = @(ww)gl_spec(ww,opt.blks);
prob.dnorm    = @(ww)gl_dnorm(ww,opt.blks);
prob.obj      = @objdalgl;
prob.softth   = @gl_softth;
prob.stopcond = opt.stopcond;
prob.ll       = min(0,yy);
prob.uu       = max(0,yy);
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
  fAT = @(x)At*x;
  fA  = @(x)A*x;
  clear At;
elseif iscell(A)
  fA  = A{1};
  fAT = A{2};
  mm = A{3};
  nn = A{4};
else
  error('A must be either numeric or cell {@(x)A*x, @(y)(A''*y), mm, nn}');
end

prob.mm       = mm;
prob.nn       = nn;

if isempty(bias)
  B = [];
else
  B = ones(mm,1);
end

[ww,bias,status]=dal(prob,ww,bias,fA,fAT,B,lambda,opt);


if all(opt.blks==opt.blks(1))
  ns=opt.blks(1);
  nc=length(ww)/ns;
  ww=reshape(ww, [ns,nc]);
end
