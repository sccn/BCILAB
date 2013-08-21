function exp_declare_symbols(varargin)
% Conveniently declare one or more symbols, to construct expressions.
%
% Example:
%    exp_declare_symbols sin cos plus pi
%    sin(cos(plus(3,pi()))) --> gives an expression that can be evaluated
%
%    the above declaration is equivalent to writing:
%    sin = exp_symbol('sin'); cos = exp_symbol('cos'); pi = exp_symbol('pi');
%
%    exp_declare_symbols f g h
%    f(g(h,h)) --> gives an expression that can only be evaluated after
%                  substituting/assigning f,g,h with other expressions
%
% See also:
%   exp_symbol
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-19

if ~iscellstr(varargin)
    error('Supplied arguments must be names for symbols.'); end

for i=1:length(varargin)
    assignin('caller',varargin{i},exp_symbol(varargin{i})); end