function str = exp_fullform(v)
% Get an (usually executable) string representation of an expression.
% String = exp_fullform(Expression)
%
% In:
%   Expression : an expression (data structure)
%
% Out:
%   String : string form of the expression
%
% Examples:
%   % generate some processed data set, then display the command that reproduces it
%   raw = flt_resample(io_loadset('data:/test.mat'))
%   tmp = flt_ica(raw,'infomax')
%   tmp = flt_iir(tmp,[0.5 2],'highpass')
%   exp_fullform(tmp)
%
%   % as previous line, but using a short-hand notation
%   ef tmp
%
% See also:
%   ef, hlp_tostring
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-15

str = hlp_tostring(v);