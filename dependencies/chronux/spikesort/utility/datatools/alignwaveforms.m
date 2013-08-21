function [offsets, aligned] = ...
	                 alignwaveforms(template, waveforms, mx)
%ALIGNWAVEFORMS    Finds best alignment of waveforms to a reference.
%   OFFSETS = ALIGNWAVEFORMS(TEMPLATE,WAVEFORMS,MAXSHIFT), where TEMPLATE
%   is an 1 x T vector and WAVEFORMS is an N x T matrix, returns an N x 1
%   vector OFFSETS giving a shift in the range [-MAXSHIFT,MAXSHIFT] for
%   each row of WAVEFORMS.  If OFFSETS(k) is X, linear interpolation of
%   WAVEFORMS(k,:) on an abcissa shifted by X samples will minimize the
%   mean squared error between WAVEFORMS(k,:) and TEMPLATE.  
%
%   If TEMPLATE is an 1 x T x D array and WAVEFORMS is an N x T x D array,
%   ALIGNWAVEFORMS finds OFFSETS that minimize the mean squared error of
%   simultaneously shifting all pages of WAVEFORMS for a given row by the
%   same amount. 
%
%   OFFSETS = ALIGNWAVEFORMS(TEMPLATE,WAVEFORMS) uses MAXSHIFT = 1.
%
%   [OFFSETS,ALIGNED] = ALIGNWAVEFORMS(TEMPLATE,WAVEFORMS,...) also
%   returns an N x (T-2*MAXSHIFT) x D matrix ALIGNED in which the rows of
%   WAVEFORMS have been shifted using linear interpolation by the amounts
%   given in OFFSETS.  Note that ALIGNED has 2*MAXSHIFT fewer columns than
%   WAVEFORMS since the interpolation operation invalidates edge samples.
%
%   The algorithm used in ALIGNWAVEFORMS actually shifts the TEMPLATE
%   waveform rather than the WAVEFORMS.  This significantly speeds up the
%   estimation but can provide erroneous results if TEMPLATE and WAVEFORMS
%   do not oversample the underlying signal. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 3),  mx = 1;  end;
if (~isnumeric(mx) || mx <= 0 || mx~=round(mx))
	error('MAXSHIFT must be a positive integer.');
end

[N,T,D] = size(waveforms);
if (~isequal([1,T,D],size(template)))
	error('The second and third dimensions of TEMPLATE must match WAVEFORMS.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
testoff = (-mx:1/3:mx)';
A = vander(testoff);  A = A(:,end-2:end);
Q = inv(A' * A) * A';  % Q fits a quadratic to errors at the different offsets

valid = (mx+1:T-mx);
for off = 1:length(testoff)
	temp = permute(template, [2,3,1]);
	templateShift(off,:,:) = ...
	    permute(interp1(1:T,temp,valid+testoff(off),'linear'), [3,1,2]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Find Alignments %%%%%%%%%%%%%%%%%%%%%%%%%%
scores = zeros(length(testoff),N);   offsets = zeros(N,1);
for n = 1:N
	for off = 1:length(testoff)
		temp = templateShift(off,:,:);  test = waveforms(n,valid,:);
		scores(off,n) = sum((temp(:) - test(:)).^2);
	end
	p = Q * scores(:,n);
	offsets(n) = -p(2)/(2*p(1));
end

offsets = -offsets;   % change to waveforms shifts instead of template shifts
offsets(offsets > mx) = mx;   offsets(offsets < -mx) = -mx;  % clip to maxshift


%%%%%%%%%%%%%%%%%%%%%%%%% Realign waveforms %%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargout < 2),  return;  end;

[X,Y] = meshgrid(1:T, 1:N);
aligned = zeros(N, T-2*mx, D);
for d = 1:D
	Xoff = X(:,valid) + repmat(offsets, [1,T-2*mx]);
	aligned(:,:,d) = interp2(X,Y,waveforms(:,:,d), Xoff, Y(:,valid), 'linear');
end


%%%%% The following code demonstrates the near-but-not-quite similar
%%%%% effects of shifting a and comparing to b vs shifting b and comparing
%%%%% to a.  When the sampling rate of a and b is oversampled relative to
%%%%% Nyquist for the underlying waveform, the difference is small.
% D = 31;   offsets = [-1:0.1:1]';  U = 100;
% up = conv2(randn(1,D*U), gausskernel([0 2*U], U) , 'same');
% b = up(35:U:end);   a = up(2:U:end);
% Sab = zeros(length(offsets),1);  Sba = zeros(length(offsets),1);
% for k = 1:length(offsets)
% 	at = interp1(1:D,a,[2:D-1]+offsets(k),'linear');   bt = b(2:end-1);
% 	Sab(k) = sum([at(:)-bt(:)].^2);
% 	
% 	bt = interp1(1:D,b,[2:D-1]+offsets(k),'linear');   at = a(2:end-1);
% 	Sba(k) = sum([at(:)-bt(:)].^2);
% end
% subplot(2,1,1);  plot([a;b]');
% subplot(2,1,2);  plot(offsets, [Sab flipud(Sba)]);

%   Last Modified By: sbm on Wed Mar  1 21:18:50 2006


%   Last Modified By: sbm on Wed Mar  1 21:19:15 2006


%   Last Modified By: sbm on Wed Mar  1 21:19:46 2006

