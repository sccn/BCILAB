% dallren - DAL with logistic loss and the Elastic-net regularization
%
% Overview:
%  Solves the optimization problem:
%   [xx, bias] = argmin sum(log(1+exp(-yy.*(A*x+bias)))) + lambda*sum(theta*abs(x)+0.5*(1-theta)*x.^2)
%
% Syntax:
%  [xx,bias,status]=dallren(xx, bias,A, yy, lambda, <opt>)
%
% Inputs:
%  xx     : initial solution ([nn,1])
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
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [ww,bias,status]=dallren(ww,bias, A, yy, lambda, theta, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt,'solver','nt',...
                     'stopcond','pdg');



prob.floss    = struct('p',@loss_lrp,'d',@loss_lrd,'args',{{yy}});
prob.fspec    = @(ww)en_spec(ww, theta);
prob.dnorm    = @(ww)en_dnorm(ww, lambda, theta);
prob.obj      = @objdalen;
prob.softth   = @en_softth;
prob.stopcond = opt.stopcond;
prob.ll       = min(0,yy);
prob.uu       = max(0,yy);
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

if isempty(bias)
  B = [];
else
  B = ones(mm,1);
end

[ww,bias,status]=dal(prob,ww,bias,fA,B,lambda,opt);

