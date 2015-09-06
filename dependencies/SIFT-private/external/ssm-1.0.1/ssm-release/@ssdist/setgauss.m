function [d y] = setgauss(d, varargin)

%@SSDIST/SETGAUSS Update the Gaussian approximation.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3
    eps     = varargin{1};
else
    y       = varargin{1};
    theta   = varargin{2};
    if any(d.type)
        eps             = y - theta;
        eps(isnan(y))   = 0; % For exponential family distributions missing data only affects ytilde
    end
end

for i = 1 : length(d.type)
    diagmask                    = d.diagmask(:, i);
    mmask                       = false(size(d.ssmat));
    mmask(diagmask, diagmask)   = true;
    if d.type(i) % Additive noise
        d.ssmat     = setdvec(d.ssmat, d.matf{i}(eps(diagmask, :)), mmask);
    else % Exponential family
        [dsubvec ytilde]    = d.matf{i}(y(diagmask, :), theta(diagmask, :));
        y(diagmask, :)      = ytilde;
        d.ssmat             = setdvec(d.ssmat, dsubvec, mmask);
    end
end

