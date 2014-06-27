% P = potCat(s,[type],[z],pots,ids) - Concatenate potential functions
%
% type and z are optional
% pots is a cell array containing the potentials e.g. pots = {'pot1', 'pot2'};
% ids  is a cell array containing the indices of the potentials e.g.
%                                                     ids = {1:100,101:200};
% s (and z if present) are split using the indices during the evaluation
%                                              feval(pot1,s(ids{1}) and so on
%
%   See also POTFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 11

function P = potCat(s,varargin) % [type], [z], pots, ids

if nargin<3
  error('wrong number of input arguments')
else
  pots = varargin{end-1}; ids = varargin{end};
end
q = 0; for i=1:length(pots), q = max(q,max(ids{i})); end

if nargin==3
  if numel(s)
    type = 'VB';      z = [];
  else                                              % return the potential types
    P.ids = ids; P.pots = pots; return
  end
elseif nargin==4
  type = varargin{1}; z = [];
elseif nargin==5
  type = varargin{1}; z = varargin{2};
else
  error('wrong number of input arguments')
end

if numel(z)==0 && strcmp(type,'VB'), P = zeros(q,4); else P = zeros(q,3); end
for i=1:length(pots)
  if numel(z)==0 && strcmp(type,'VB')
    P(ids{i},:) = feval(pots{i}, s(ids{i}), type, []);
  else
    P(ids{i},:) = feval(pots{i}, s(ids{i}), type, z(ids{i}));
  end
end