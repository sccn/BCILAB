function [A Pe_out err Pa_out] = DEKF(y,p,lambda,verb,downsampleFactor,REINIT)
%Dual Extended Klman Filter (DEKF) for MVAR parameter estimation
%
% Arguments:
% A: Estimated time-varying parameters, A = [A1 A2 ... Ar]
% Pe: Process noise covariance matrix
% e:  residuals
% Pa: MVAR parameter covariance matrix
% y: (M x LEN) data matrix
% p: Model order
%
%
% References:
% [1] E. A. Wan, and A. T. Nelson, “Neural dual extended Kalman filtering: applications in speech
% enhancement and monaural blind signal separation,” in Neural Networks for Signal Processing
% of the 1997 IEEE Workshop, 1997, pp. 466-475.
%
% [2] E. A. Wan, and A. T. Nelson, "Dual Extended Kalman Filter Methods," Kalman Filtering and
% Neural Networks, pp. 123-173: John Wiley & Sons, Inc., 2002.
%
% [3] S. Haykin, Kalman Filtering and Neural Networks, p.^pp. 304: John Wiley and Sons, 2001.
%
%
% See also: 'Linear Kalman Filter' MATLAB implementation, written by Amir Omidvarnia, available at:
% http://www.mathworks.com/matlabcentral/fileexchange/29127-linear-kalman-filter
%
% Written by: Amir Omidvarnia
% Modified/Optimized by: Tim Mullen

%%

if nargin<2
    error('SIFT:DEKF:Unsufficient arguments');
end

if nargin<3
    % update coefficent for parameter covariance matrix (forgetting factor)
    lambda = .02;
end

if nargin<4
    verb = 1;
end

if nargin<5
    downsampleFactor = 1;
end

if nargin<6
    % if REINIT=true, then reinitialize params between trials.
    % Otherwise, start filtering each trial from previous solution
    REINIT = false;
end

if verb
    h=waitbar(0,sprintf('fitting VAR[%d] model [mode=%s] ...', ...
        p, 'Dual EKF'));
end

%%
[M LEN NTR] = size(y);
% M = size(y,1);                            % Number of states (here, M = N)
% LEN = size(y,2);                          % Number of the multivariate observations


dslen = ceil((LEN-downsampleFactor+(downsampleFactor>1))/downsampleFactor);

% Initialize output variables
A       = zeros(M,M*p,dslen,NTR);
if nargout>=2
    Pe_out  = zeros(M,M,dslen,NTR);                % (EKF 1) Process noise covariance matrix (TM)
end
if nargout>=3
    err = zeros(size(y,1),dslen,NTR);
end
if nargout>=4
    Pa_out  = zeros(M*M*p,M*M*p,dslen,NTR);        % (EKF 1) Parameter covariance matrix (TM)
end     

for tr = 1 : NTR
    
    %% Initial parameters for Dual Extended Kalman Filter
    
    if tr==1 || REINIT
        
        
        %%%%% (EKF 1)
        xh = zeros(M*p,LEN);                      % (EKF 1) Initial a-posteriori states (Mp x 1)
        Px = .1*eye(M*p);                         % (EKF 1) Initial a-posteriori state covariance matrix
        R = eye(M);
        B = zeros(M*p,M);                        % (EKF 1) Relationship between states and the process noise ( x(k) = F[x(k-1)] + B*v(k) )
        B(1:M,:) = eye(M);                       % (EKF 1) B = [I 0 ... 0]'
        
        %%%%% EKF 2
        Pa = eye(M*M*p);                         % (EKF 2) Initial a-posteriori parameters covariance matrix
        % Ah_ = sparse([],[],[],M*p,M*p,M*M*p+M*(p-1));   % (EKF 2) Initial a-posteriori parameters estimates (matrix form of 'ah' plus identity matrices)
        Ah_ = zeros(M*p,M*p);
        Ah_(1:M,1:M*p) = .1*randn(M,M*p);       % (EKF 2) Initial a-posteriori parameters estimates (matrix form of 'ah' plus identity matrices)
        Ah_(M+1:M*p,1:(p-1)*M) = eye(M*(p-1)); %speye(M*(p-1));
        
        % for r = 2 : p
        %     Ah_((r-1)*M+1:r*M,(r-2)*M+1:(r-1)*M) = eye(M);
        % end
        
        for r = 1 : p
            xh((r-1)*M+1:r*M,p+1) = y(:,p-r+1,tr);
        end
        
        %%%%% Mutual variables between EKF 1 and DEKF 2
        Q = 10*eye(M);                            % (EKF 1,2) Initial process noise covariance matrix
        C = B';                                   % (EKF 1,2) Measurement matrix (identity matrix, C = B')
        
        I_CHp       = eye(M*p);
        I_CH2p      = eye(M*M*p);
        % I_CHp_m1    = eye(M*(p-1));
        
        %% DEKF starts ....
        % Ah_ = Ah;                                 % Ah_(k) = Ah(k-1)
        
    end
   

    % store default covariance matrix estimates for first p samples
    Pe_out(:,:,1:ceil(p/downsampleFactor),tr) = repmat(Q,[1,1,ceil(p/downsampleFactor)]);
    
    curval = ceil(p/downsampleFactor)+1;
        
    for n = p+1 : LEN
        
        if verb==2 && ~mod(n,floor(LEN/10))
            waitbar(n/LEN,h,...
                {sprintf('Trial (%d/%d)',tr,NTR), ...
                sprintf('fitting VAR[%d] model [mode=%s] (%d/%d) ...',p,'Dual EKF',n,LEN)});
        end
        
        
        
        if any(isnan(y(:,n,tr)))
            % don't update on missing data
            continue;
        end
        
        [J_x J_A] = MVAR_JacCSD(Ah_,xh(:,n-1),p);  % xh(k) = F(A(k-1) * xh(k-1)) = Ah(k-1) * xh(k-1)
        
        %% EKF 1 ---> States estimation
        %---------- Time Update (EKF1) ----------
        Rv = B * Q * B';                          % According to Haykin's book
        xh_ = Ah_ * xh(:,n-1);                    % xh_(k) = A_h(k-1) * xh(k-1)
        Px_ = J_x * Px * J_x' + Rv;               % Px_(k) = A_h(k-1) * Px(k-1) * A_h(k-1)' + B * Q * B'
        
        %---------- Measurement Update (EKF1) ----------
        Rn = R; %R(:,:,n-1);                      % According to Haykin's book
        Kx = Px_ * C' /(C * Px_ * C' + Rn);       % Kx(k)  = Px_(k) * C' * inv(C * Px_(k) * C' + R)
        Px = (I_CHp - Kx * C) * Px_;          % Px(k)  = (I - Kx(k) * C) * Px_(k)
        e = y(:,n,tr) - C * xh_;                     % inov(k) = y(k) - C * Ah_(k) * xh(k-1)
        xh(:,n) = xh_ + Kx * e;                   % xh(k)  = xh_(k) + Kx(k) * (y(k) - C * xh_(k))
        
        %% EKF 2 ---> Parameters estimation
        %---------- Time Update (EKF2) ----------
        ah_ = reshape(Ah_(1:M,:)',M*M*p,1);    % ah_ = vec(Ah(k-1))
        Rr = lambda*Pa;                           % Rr = lambda * Pa(k-1)
        Pa_ = Pa + Rr;                            % Pa_(k) = Pa(k-1) + Rr
        
        %---------- Measurement Update (EKF2) ----------
        %%%%%% Compute DfDa and H (H(k) = C*DfDa(k-1))
        H = C * (-J_A);                           % J_A = -DfDa;
        
        %%%%%%%%%%%%%%%
        Re = (Rn + Q);                            % According to Haykin's book
        Ka = Pa_ * H' /(H * Pa_ * H' + Re);       % Ka(k) = Pa_(k) * H(k) * inv(H(k) * Pa_(k) * H(k)' +R + Q)
        
        %% NOTE: This line is takes the most time -- optimize if possible (GPU?)
        Pa = (I_CH2p - Ka * H) * Pa_;           % Pa(k) = (I - Ka(k) * H(k)) * Pa_(k)
        %% END NOTE
        
        ah = ah_ + Ka * (y(:,n,tr)- C * Ah_ * xh(:,n-1));    % ah(k) = ah_(k) + Ka(k) * (y(k) - yh_(k))
        
        %---------- Re-arrange vector ah(k) into the matrix Ah(k) ----------
        Ah_(1:M,:) = reshape(ah,M*p,M)';
        
        %     Ah(M+1:M*p,1:(p-1)*M) = I_CHp_m1;
        %     for r = 2 : p
        %         Ah((r-1)*M+1:r*M,(r-2)*M+1:(r-1)*M, n) = eye(M);
        %     end
        
        if ~mod(n-1,downsampleFactor)
            A(:,:,curval,tr)        = Ah_(1:M,:); % Estimated time-varying parameters, A = [A1 A2 ... Ar]
            if nargout>=2
                Pe_out(:,:,curval,tr)   = (1-lambda)*Q + lambda*(e*e');            % noise covariance matrix
            end
            if nargout>3
                err(:,n,tr) = e;
            end
            if nargout>=4
                Pa_out(:,:,curval,tr)   = Pa;
            end
            curval = curval+1;
        end
    end
    
end

if verb, close(h); end
