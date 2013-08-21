function [f,F] = pulseTrain( xx, a, b, c, d, r, S, multipleOut )
% f = pulseTrain( x, a, b, c, d )
%   returns the value of a trapezoidal function evaluated at x
%   (as a row vector)
% f = pulseTrain( x, a, b, c, d, R )
%   assumes the trapezoidal pulse repeats every "R" seconds
% f = pulseTrain( x, a, b, c, d, R, S)
%   scales the output by S (otherwise, the pulse has unit height)
% [f,F] = pulseTrain( ... )
%   returns a matrix F, where each row corresponds to the pulse
%   from a different repetition
%
% if b < a, then assumes the input is in the form:
%   a = a (offset)
%   b = a + b ( rise time )
%   c = c + b ( pulse duration )
%   d = d + c ( fall time )
% (or, if b is negative, it will force this mode )
%
%   The window looks like:
%            "b"          "c"
%             * --------- *
%           *              * 
%   ______ *                 * ________
%         "a"                "d"
%   <--------------- r --------------->

if nargin > 5 && ~isempty(r)
    x = rem(xx,r);
end
if nargin < 7 || isempty(S), S = 1; end

if b < a || b < 0
    b = abs(b);
    b  = a + b;
    c = c + b;
    d = d + c;
end

if numel(x) == 1
    % a scalar input, so calculation is very simple:
    if x < a || x > d, f = 0; return; end
    if x > b && x < c, f = S*1; return; end

    if x > c % and we know x < d
        x = x - c;
        f = S*( 1 - x/(d-c) );
    else
        % we know a < x < b
        x = x-a;
        f = S*( x/(b-a) );
    end
else
    
    % vectorize it
    ba = b-a; dc = d-c;
    
    f = (x<d).*(x>a).*( ...
            (x<=b).*(x-a)/ba + ...
            (x>b).*( ...
                (x<=c) + ...
                (x>c).*(1-( x-c )/dc) )   );
    f = S*f;
        
end

if nargout >= 2
    T = xx(end);
    n = length(x);
    nRep = ceil( T / r );
    F = zeros( nRep, n);
    for n = 1:nRep
        F(n,:) = f.*( ( xx >= (n-1)*r ) & ( xx < n*r ) );
    end
end
    
    