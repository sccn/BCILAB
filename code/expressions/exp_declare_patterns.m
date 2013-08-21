function exp_declare_patterns(varargin)
% Conveniently declare one or more pattern expressions, for use in pattern matching.
% exp_declare_patterns x y z
%
%  the above is equivalent to writing:
%  x = exp_pattern(exp_symbol('x')); y = exp_pattern(exp_symbol('y')); z = exp_pattern(exp_symbol('z'));
%
% See also:
%   exp_pattern
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-19

if ~iscellstr(varargin)
    error('supplied arguments must be names for symbols.'); end

for i=1:length(varargin)
    assignin('caller',varargin{i},exp_pattern(exp_symbol(varargin{i}))); end