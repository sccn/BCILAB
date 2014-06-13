% Penalised least squares solved by CONJUGATE GRADIENTS using the function
% minimize.m [1] by Carl E. Rasmussen.
%
%   [u,phi] = plsCG(u,X,y,B,t,opt,lam,pen,varargin)
%  
%   See also PLSSOLVERS.M.
%
% (c) by Hannes Nickisch, Philips Research, 2013 August 30
%
% [1] http://www.kyb.tuebingen.mpg.de/bs/people/carl/code/minimize/

function [u,phi] = plsCG(u,X,y,B,t,opt,lam,pen,varargin)

% optimisation parameters
if isfield(opt,'nMVM')                                 % maximal number of steps
  nMVM = opt.nMVM;
else
  nMVM = 100;
end

% information parameters
if isfield(opt,'output')              % flag saying whether some output is shown
  output = opt.output;
else
  output = false;
end

OCT = exist('OCTAVE_VERSION') ~= 0;           % check if we run Matlab or Octave
if output || OCT                           % evalc is not implemented for Octave
  [u,phi] = minimize(u,'phi',-nMVM,X,y,B,t,lam,pen,varargin{:});
else
  [txt,u,phi] = evalc('minimize(u,''phi'',-nMVM,X,y,B,t,lam,pen,varargin{:})');
end
phi = phi(end);