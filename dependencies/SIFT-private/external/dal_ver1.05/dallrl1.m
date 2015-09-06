% dallrl1 - DAL with logistic loss and the L1 regularization
%
% Overview:
%  Solves the optimization problem:
%   [xx, bias] = argmin sum(log(1+exp(-yy.*(A*x+bias)))) + lambda*||x||_1
%
% Syntax:
%  [ww,bias,status]=dallrl1(ww0, bias0, A, yy, lambda, <opt>)
%
% Inputs:
%  ww0    : initial solution ([nn,1])
%  bias0  : initial bias (set [] if bias term is unnecessary)
%  A      : the design matrix A ([mm,nn])
%  yy     : the target label vector (-1 or +1) ([mm,1])
%  lambda : the regularization constant
%  <opt>  : list of 'fieldname1', value1, 'filedname2', value2, ...
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
% m = 1024; n = 4096; k = round(0.04*n); A=randn(m,n);
% w0=randsparse(n,k); yy=sign(A*w0+0.01*randn(m,1));
% lambda=0.1*max(abs(A'*yy));
% [ww,bias,stat]=dallrl1(zeros(n,1), 0, A, yy, lambda);
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [ww,bias,status]=dallrl1(ww,bias, A, yy, lambda, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt,'solver','cg',...
                     'stopcond','pdg');

if ~isequal(unique_bc(yy), [-1;1])
  error('yy must be a column vector of -1''s and 1''s');
end



prob.floss    = struct('p',@loss_lrp,'d',@loss_lrd,'args',{{yy}});
prob.fspec    = @(xx)abs(xx);
prob.dnorm    = @(vv)max(abs(vv));
prob.obj      = @objdall1;
prob.softth   = @l1_softth;
prob.stopcond = opt.stopcond;
prob.ll       = min(0,yy);
prob.uu       = max(0,yy);
prob.Ac       =[];
prob.bc       =[];
prob.info     =[];

if isequal(opt.solver,'cg')
  prob.hessMult = @hessMultdall1;
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

if isempty(bias)
  B = [];
else
  B = ones(mm,1);
end


[ww,bias,status]=dal(prob,ww,bias,fA,B,lambda,opt);

