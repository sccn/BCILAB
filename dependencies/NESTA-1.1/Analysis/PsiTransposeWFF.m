function y = PsiTransposeWFF(x,w_type,log_length,min_scale,max_scale,shift_redundancy,freq_redundancy, PLOT)
% Windowed Fourier Frame Analysis
% y = PsiTransposeWFF( x, w_type, log_length, min_scale, max_scale,...
%       shift_redundancy, freq_redundancy, plot )
%   w_type is the type of window.  Currently, this supports 'isine' (iterate sine)
%   and 'gaussian'.  Default: 'isine' (to use default, set this to [] ).
% If plot == true, then will make intermediate plots (default: false)
%   to use this feature, set the global variable "dT" to the Nyquist time
%   (note: units are probably not correct).
% Core code written by Peter Stobbe; modifications by Stephen Becker
if isempty(w_type), w_type = 'isine'; end
if nargin < 8 || isempty(PLOT), PLOT = false; end

% w is a is a vector of with the window of the largest scale
% smaller scale windows are just this subsampled
[w, window_redundancy] = make_window(max_scale,w_type);

y = [];
c = ((max_scale - min_scale + 1)*window_redundancy*2.^((1:max_scale)'+freq_redundancy+shift_redundancy)).^-.5;

% % If parallel computing toolbox exists, do this in parallel:
% if exist('parfor','builtin')
%   parfor k = min_scale:max_scale
%     M = 2^(log_length-k) +(2^(log_length-k)+1)*(2^shift_redundancy-1);
%     z = [myRepMat(x,2^shift_redundancy); zeros(2^k - 2^(k-shift_redundancy),2^shift_redundancy)];
%     z = reshape(z,2^k,M);
%     z = z.*myRepMat(w(2^(max_scale-k)*(1:2^k)'),M);
%     z(2^(k+freq_redundancy),M) = 0;
%     z = fft(z);
%     z = [z(1,:)*c(k);       real(z(2:end/2,:))*c(k-1); ...
%         z(end/2+1,:)*c(k); imag(z(end/2+2:end,:))*c(k-1)];
%     y = [y; z(:)];
%   end
% else
    
    
    for k = min_scale:max_scale
        M = 2^(log_length-k) +(2^(log_length-k)+1)*(2^shift_redundancy-1);
        z = [myRepMat(x,2^shift_redundancy); zeros(2^k - 2^(k-shift_redundancy),2^shift_redundancy)];
        %     z = [repmat(x,1,2^shift_redundancy); zeros(2^k - 2^(k-shift_redundancy),2^shift_redundancy)];
        z = reshape(z,2^k,M);
        z = z.*myRepMat(w(2^(max_scale-k)*(1:2^k)'),M);
        %     z = z.*repmat(w(2^(max_scale-k)*(1:2^k)'),1,M);
        z(2^(k+freq_redundancy),M) = 0;
        z = fft(z);
        z = [z(1,:)*c(k);       real(z(2:end/2,:))*c(k-1); ...
            z(end/2+1,:)*c(k); imag(z(end/2+2:end,:))*c(k-1)];
        % Stephen is adding this plotting feature:
        if PLOT
            global dT
            if isempty(dT), disp('global variable dT must be set!'); dT = 1; end
            figure(k-min_scale+1); clf
            colormap(1-gray(256))
            fprintf('Size of transformation at level %d is %d x %d\n',...
                k-min_scale+1,size(z,1), size(z,2) );
            yy = size(z,1);
            xx = size(z,2);
            zz = z(1:yy/2, :);
            mx = max(max(abs(zz)));
            image( linspace(0,dT*2^log_length,xx), linspace(0,1/dT/2,yy/2), 256*abs(zz)/mx );
            %xlim([0,2^log_length * dT * shift_redundancy]);
            xlabel('time in s'); ylabel('frequency in Hz');
            title('Time Frequency plot.  Warning: units may not be correct')
            colorbar;
        end
        y = [y; z(:)];
    end
% end
%---------------


% Stephen is adding this on 8/19/08:
% for small n, the overhead in repmat is actually
% significant, since the fft is so fast already
function B = myRepMat(A,n)
B = A(:,ones(n,1));

