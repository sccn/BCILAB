function [theta ksivar] = htd(d, D, s, Phi, Theta, etavar, tol)

%HTD Hillmer-Tiao decomposition.
%   [theta ksivar] = HTD(d, D, s, Phi, Theta[, etavar, tol])
%       d is the order of ordinary differencing.
%       D is the order of seasonal differencing.
%       s is the seasonal period.
%       Phi contains the autoregressive factor for each component.
%       Theta is the moving average factor.
%       etavar is the disturbance variance.
%       tol is the tolerance for rounding error (which causes imaginary
%           numbers or negative covariances).
%       theta is the decomposed moving average factor for each component.
%       ksivar is the decomposed disturbance variance for each component
%           (including irregular).

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 7, tol = 10^3*eps; end
if nargin < 6, etavar = 1; end
if isempty(Phi), Phi = {}; elseif isnumeric(Phi), Phi = {Phi}; end

eps_resolution  = 0.0001;

% Generate differencing and sum factors
DPoly   = [-1 1];
for i = 2 : d, DPoly = [0 DPoly] + [-DPoly 0]; end
UPoly   = ones(1, s);
for i = 2 : D, UPoly = conv(UPoly, ones(1, s)); end

% Construct total and individual autoregression factors
% Number of autoregression factors = two differencing factors + AR factors
m       = length(Phi);
phi     = Phi;
if D > 0
    m       = m + 1;
    phi     = {UPoly(end-1:-1:1) phi{:}};
end
if d > 0
    m       = m + 1;
    phi     = {DPoly(end-1:-1:1) phi{:}};
end
p       = zeros(1, m);
for i = 1:m
    phi{i}  = [phi{i}(end:-1:1) 1];
    p(i)    = length(phi{i});
end
Theta       = [Theta(end:-1:1) 1];
q           = length(Theta);

n           = max([p q]);
S           = zeros(n);
S(n, n)     = 2;
S(n-1, n-1) = 1;
for i = n-2:-1:1, S(i, :) = [S(i+1, 2:n) 0] - S(i+2, :); end
S(n, n)     = 1;

remphi      = cell(1, m);
for i = m : -1 : 1
    phi_sym = conv(phi{i}, phi{i}(p(i):-1:1));
    phi{i}  = phi_sym(1:p(i))*S(n-p(i)+1:n, n-p(i)+1:n);
    if i == m, remphi{i} = phi{i};
    else remphi{i} = conv(remphi{i+1}, phi{i});
    end
end
Theta_sym   = conv(Theta, Theta(q:-1:1));
Theta       = Theta_sym(1:q)*S(n-q+1:n, n-q+1:n);

% Polynomial division
[Q_eps R]   = polydiv(Theta, remphi{1});

% Partial fraction
Q           = cell(1, m);
for i = 1 : m-2
    [Q{i} R]    = polypartfrac(phi{i}, remphi{i+1}, R);
end
[Q{m-1} Q{m}]   = polypartfrac(phi{m-1}, phi{m}, R);
Q{m}            = [0 Q{m}]; % Pad 0 to make Q{m} the same 'degree' as phi{m}

% Add extra pure MA component
if ~isscalar(Q_eps)
    m       = m + 1;
    p(m)    = length(Q_eps);
    phi{m}  = [zeros(1, length(Q_eps)-1) 1];
    Q{m}    = Q_eps;
end

% Canonical decomposition
mineps  = zeros(1, m);
X       = 2*cos(0:eps_resolution*pi:pi);
warn_S  = warning('off', 'MATLAB:divideByZero'); % The spectrum will normally have spikes
for i = 1:m, mineps(i) = min(polyval(Q{i}, X)./abs(polyval(phi{i}, X))); end
warning(warn_S.state, 'MATLAB:divideByZero');

% The minimums of the spectrums of each part could be invalid for numerous
% reasons:
% 1. -Inf most likely occurs in trend and seasonal differencing parts,
%    where the denominator (phi) have unit roots, if the numerator are
%    negative at those frequencies no admissable decomposition is possible.
% 2. NaN occurs for 0/0 or +/-Inf/Inf, it may be theoretically possible to
%    determine whether the result is 0 or Inf, but since MATLAB ignores NaN
%    when taking the minimum this means all frequencies are NaN, hence
%    computationally 'inadmissible'.
% 3. Inf most likely occurs when one of the partial fraction has no
%    solution, the quotient and remainder then becomes infinite
%    polynomials.
invalid_mineps  = (mineps == -Inf | isnan(mineps) | mineps == Inf);

% Factorization
theta   = cell(1, m);
ksivar  = zeros(1, m+1);
for i = 1 : m
    thetafactor = Q{i} - mineps(i)*phi{i};
    if thetafactor(1) == 0 || invalid_mineps(i)
        %%%% TODO: This is a temporary solution that outputs a 'numerical'
        %%%% result without generating an error. Not all functions take NaN
        %%%% so this just postpones the inevitable.
        theta{i}    = repmat(NaN, 1, length(phi{i})-1);
        ksivar(i)   = 0;
    else
        r           = polyrootpairs(thetafactor, tol);
        if any(isnan(r))
            %%%% TODO: See above
            theta{i}    = repmat(NaN, 1, length(phi{i})-1);
            ksivar(i)   = 0;
        else
            theta{i}    = [zeros(1, nnz(r==Inf)) poly(r)];
            theta{i}    = theta{i}/theta{i}(end);
            theta_sym   = conv(theta{i}, theta{i}(end:-1:1));
            ksivar(i)   = thetafactor(end)/(theta_sym(1:p(i))*S(n-p(i)+1:n, n));
            theta{i}    = theta{i}(end-1:-1:1);
        end
    end
end

% Irregular variance
if isscalar(Q_eps)
    ksivar(m+1) = Q_eps + sum(mineps);
else
    ksivar(m+1) = sum(mineps);
end
if ksivar(m+1) < 0 || isnan(ksivar(m+1))
    warning('ssm:htd:inadmissible', 'No acceptable decomposition exists.');
    ksivar(m+1) = NaN;
end
ksivar      = etavar*ksivar;
    
function [q r] = polydiv(a, b)
% a(x) = q(x)*b(x) + r(x)
na  = length(a); % >= 1
nb  = length(b); % >= 1
if na >= nb && nb > 1
    nq  = na - nb + 1;
    q   = zeros(1, nq);
    for i = 1 : nq
        q(i)        = a(i)/b(1);
        a(i:i+nb-1) = a(i:i+nb-1) - q(i)*b;
    end
    r   = a(nq+1:na);
elseif nb == 1
    q   = a/b;
    r   = 0;
else % nb > 1 && na < nb
    q   = 0;
    r   = a;
end

function [s t] = polypartfrac(a, b, c)
% s(x)/a(x) + t(x)/b(x) = c(x)/a(x)b(x)
% where a(x) has degree p, b(x) has degree q, and c(x) should have degree p+q-1
p   = length(a) - 1;
q   = length(b) - 1;
A   = zeros(p+q);
for i = 1 : p+q
    if i <= p, A(:, i) = [zeros(i-1, 1); b'; zeros(p-i, 1)];
    else j = i - p; A(:, i) = [zeros(j-1, 1); a'; zeros(q-j, 1)];
    end
end
if length(c) < p+q, c = [zeros(1, p+q-length(c)) c]; end
st  = inv(A)*c';
s   = [0 st(1:p)']; % Pad 0 to make s the same 'degree' as a
t   = st(p+1:p+q)';

function r = polyrootpairs(d, tol)
% r is the roots of c on or outside the unit circle

% Find inverse pair roots
n       = length(d);
A       = diag(ones(1, n-2), -1);
A(1, :) = -d(2:n)/d(1);
lambda  = eig(A);
if isreal(lambda)
    r   = (lambda + sqrt(lambda.^2 - 4))/2; % roots always inverse pairs, so only need one of them
else
    r           = zeros(size(lambda));
    relambda    = (imag(lambda) == 0)';
    for i = find(~relambda)
        a       = real(lambda(i));
        b       = imag(lambda(i));
        k       = a^2 - b^2 - 4;
        m       = 2*a*b;
        h       = realsqrt((abs(k) + realsqrt(k^2 + m^2))/2);
        if k >= 0
            c   = h;
            d   = m/(2*c);
        else
            d   = sign(m)*h;
            c   = m/(2*d);
        end
        r(i)    = complex(a + c, b + d)/2;
    end
    r(relambda) = (lambda(relambda) + sqrt(lambda(relambda).^2 - 4))/2; % roots always inverse pairs, so only need one of them
end
r(abs(imag(r))<realsqrt(tol))   = real(r(abs(imag(r))<realsqrt(tol))); % Get rid of imaginary parts from rounding error
absr        = abs(r);
r(absr<1)   = 1./r(absr<1); % inverse some so that all will have modulus > 1
ucmask      = (abs(1-absr)<tol) & (imag(r) ~= 0);
if any(ucmask)
    % Make sure unit roots that are complex are paired correctly
    % (there will always be an even number of them)

    %%%% TODO: There's actually sometimes an odd number of unit roots, due
    %%%%       to some of lambda being real and absolute value < 2, say
    %%%%       2 - 10^-8 or something, this results in a single complex
    %%%%       root with no hope of having a conjugate, as 10^-8 ~= eps^0.5
    %%%%       this may be due to rounding error, yet lambda comes directly
    %%%%       out of eig in line 161 ...
    if mod(nnz(ucmask), 2)
        warning('ssm:htd:RoundingError', 'Unrecoverable rounding error.');
        r   = repmat(NaN, n-1, 1);
        return;
    end

    % Strip existing conjugate pairs
    for i = find(ucmask')
        for j = i + find(ucmask(i+1:end)')
            if ucmask(i) && ucmask(j) && isreal(r(i) + r(j))
                ucmask([i j]) = false;
            end
        end
    end
    % if any(ucmask) at this point they are not conjugates
    % (length of uc will still be even)
    uc              = r(ucmask);
    uc              = [uc; conj(uc)];
    [ucind{1:2}]    = sort(angle(uc));
    uc              = uc(ucind{2});
    r(ucmask)       = [uc(1:end/4) uc(end:-1:(3*end/4)+1)]';
end
