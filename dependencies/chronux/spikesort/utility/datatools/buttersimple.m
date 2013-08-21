function [B, A] = buttersimple(Wn, option)
%BUTTERSIMPLE      Fixed order Butterworth digital filter design.
%   [B, A] = BUTTERSIMPLE(Wn), where 0.0 < Wn < 1.0, computes the
%   numerator/denominator coefficients for a 6th order lowpass digital
%   Butterworth filter.
%
%   [B, A] = BUTTERSIMPLE(Wn, 'high'), where 0.0 < Wn < 1.0, computes the
%   coefficients for a 6th order high-pass digital Butterworth filter.
%
%   [B, A] = BUTTERSIMPLE([W1 W2]), with 0.0 < W1 <1.0 and 0.0 < W2 < 1.0,
%   computes the coefficients for a 6th order bandpass digital filter with
%   passband W1 < W < W2.
%
%   Based on TMW Signal Processing Toolbox function BUTTER, but without
%   any of the generality and does not require any toolboxes.

% Step 1: Compute analog frequencies
if (~all(Wn>0 & Wn<1.0)),  
	error('Cutoff frequencies must be in the interval (0,1).');  
end
Wn = 4*tan(pi*Wn/2);
if (length(Wn) == 2)
	mode = 'B';
	Bw = diff(Wn);          % bandwidth
	Wn = sqrt(prod(Wn));	% center frequency
elseif (length(Wn) == 1)
	if (nargin > 1)
        if (~strcmpi(option, 'high')),  error('Unknown option.');   end;
        mode = 'H';
    else  mode = 'L';
	end
else
    error('Wn can not have more than two elements.');
end

% Step 2: Make Butterworth lowpass prototype
if (mode == 'B')  % if bandpass, start w/ 3rd order prototype
	a = [-1 0 0; 1 -1 -1; 0 1 0];
	b = [1 0 0]';    
    c = [0 0 1];
	d = 0;
else  % otherwise 6th order lowpass prototype
	a = diag([1 1 1 1 1], -1) + diag([-1 0 -1 0 -1],1);
    for j = [1,3,5]
        a(j,j) = 2*cos((1-(j/12))*pi);
    end
	b = [1 0 0 0 0 0]';
    c = [0 0 0 0 0 1];
	d = 0;
end
		
% Step 3: Convert to lowpass/bandpass/highpass & shift frequency
switch(mode)
    case 'L',
        a = Wn*a;
        b = Wn*b;
%       c = c;
%       d = d;
    case 'B',
        a = Wn*[a/(Wn/Bw) eye(3); -eye(3) zeros(3)];
        b = [b*Bw; 0; 0; 0];
        c = [c 0 0 0];
        d = 0;
    case 'H',
        d = d - c/a*b;        
        c = c/a;
        b = -Wn*(a\b);       
        a =  Wn*inv(a);
end

% Step 4: Find discrete equivalent
t1 = eye(6) + a/4;    t2 = eye(6) - a/4;
d = c/t2*b/4 + d;
c = c/t2/sqrt(2);
b = t2\b/sqrt(2);
a = t2\t1;

% Step 5: Convert to polynomial form
A = poly(a);
B = poly(a-b*c)+(d-1)*A;
