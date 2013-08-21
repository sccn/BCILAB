function z=hatmatrix(varargin)

fit = locfit(varargin{:},'module','hatm');
z = lfknots(fit)';

return;
