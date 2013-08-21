function x = PsiWFF(y,w_type,log_length,min_scale,max_scale,shift_redundancy,freq_redundancy);
% Windowed Fourier Frame Synthesis
% x = PsiWFF(y,w_type,log_length,min_scale,max_scale,shift_redundancy,freq_redundancy);
%   w_type is the type of window.  Currently, this supports 'isine' (iterate sine)
%   and 'gaussian'.  Default: 'isine' (to use default, set this to [] ).
% written by Peter Stobbe, modified by Stephen Becker
if isempty(w_type), w_type = 'isine'; end

% w is a is a vector of with the window of the largest scale
% smaller scale windows are just this subsampled
[w, window_redundancy] = make_window(max_scale,w_type);

x = zeros(2^log_length,1);
c = ((max_scale - min_scale + 1)*window_redundancy*2.^((1:max_scale)'+freq_redundancy+shift_redundancy)).^-.5;
cur = 0;

% % If parallel computing toolbox exists, do this in parallel:
% if exist('parfor','builtin')
%     parfor k = min_scale:max_scale
%     M = 2^(log_length-k) +(2^(log_length-k)+1)*(2^shift_redundancy-1);
%     N = 2^(k+freq_redundancy);
%     
%     % avoid the "cur" problem:
%     kk = min_scale:k;
%     Ms = 2.^(log_length-kk) +(2.^(log_length-kk)+1).*(2.^shift_redundancy-1);
%     Ns = 2.^(kk+freq_redundancy);
%     cur = sum( Ms(1:end-1).*Ns(1:end-1) );
%     
%     z = reshape(y(cur + (1:N*M)),N,M);
% %     cur = cur + N*M;
%     z = [z(1,:)*c(k);     z(2:N/2,:)*c(k-1);...
%         z(N/2+1,:)*c(k); -1i*z(N/2+2:N,:)*c(k-1)];
%     z = real(fft(z));
%     z = z(1:2^k,:).*myRepMat(w(2^(max_scale-k)*(1:2^k)),M);
%     z = reshape(z,[],2^shift_redundancy);
%     x = x + sum(z(1:2^log_length,:),2);
%     end
%     
% else
    
    for k = min_scale:max_scale
        M = 2^(log_length-k) +(2^(log_length-k)+1)*(2^shift_redundancy-1);
        N = 2^(k+freq_redundancy);
        z = reshape(y(cur + (1:N*M)),N,M);
        cur = cur + N*M;
        z = [z(1,:)*c(k);     z(2:N/2,:)*c(k-1);...
            z(N/2+1,:)*c(k); -1i*z(N/2+2:N,:)*c(k-1)];
        z = real(fft(z));
        z = z(1:2^k,:).*myRepMat(w(2^(max_scale-k)*(1:2^k)),M);
        %     z = z(1:2^k,:).*repmat(w(2^(max_scale-k)*(1:2^k)),1,M);
        z = reshape(z,[],2^shift_redundancy);
        x = x + sum(z(1:2^log_length,:),2);
    end
% end

%---------------


% Stephen is adding this on 8/19/08:
% for small n, the overhead in repmat is actually
% significant, since the fft is so fast already
function B = myRepMat(A,n)
B = A(:,ones(n,1));

