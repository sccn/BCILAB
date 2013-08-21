% Penalised least squares solved by the LIMITED memory
% BROYDEN-FLETCHER-GOLDFARB-SHANNO (L-BFGS) algorithm which is a quasi Newton or
% variable metric method. 
% We use Peter Carbonetto's "Matlab interface for L-BFGS-B" [1].
%
%   [u,phi] = plsTN(u,X,y,B,opt,lam,pen,varargin)
%
% Additional options:
%  opt.
%      nonneg:    enforce positive u                            [default  false]
%
%   See also PLSSOLVERS.M.
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 23
%
% [1] http://www.cs.ubc.ca/~pcarbo/lbfgsb-for-matlab.html

function [u,phi] = plsLBFGS(u,X,y,B,opt,lam,pen,varargin)

% optimisation parameters
if isfield(opt,'nMVM')                                 % maximal number of steps
  nMVM = opt.nMVM;
else
  nMVM = 100;
end
if isfield(opt,'LBFGSnonneg')                           % nonnegativity desired?
  nonneg = opt.LBFGSnonneg;
else
  nonneg = false;
end

% information parameters
if isfield(opt,'output')              % flag saying whether some output is shown
  output = opt.output;
else
  output = false;
end

OCT = exist('OCTAVE_VERSION') ~= 0;           % check if we run Matlab or Octave
if output || OCT                           % evalc is not implemented for Octave
  if output, if exist('fflush','builtin'), fflush(stdout); end, end
  if nonneg
    if output, fprintf('    => non-negativity enforced\n'), end
    [u,phi] = minimize_lbfgs(u,'phi', abs(nMVM),X,y,B,lam,pen,varargin{:});
  else
    [u,phi] = minimize_lbfgs(u,'phi',-abs(nMVM),X,y,B,lam,pen,varargin{:});
  end
  if output, fprintf('\n'); if exist('fflush','builtin') fflush(stdout); end, end
else
  [t,u,phi] = evalc('minimize_lbfgs(u,''phi'',-nMVM,X,y,B,lam,pen,varargin{:})');
end

phi = phi(end);