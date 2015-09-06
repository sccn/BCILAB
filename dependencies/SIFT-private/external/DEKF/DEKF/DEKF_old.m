function [A Pe e Pa] = DEKF(y,p,lambda,verb)
%Dual Extended Klman Filter (DEKF) for MVAR parameter estimation
%
% Arguments:
% A: Estimated time-varying parameters, A = [A1 A2 ... Ar]
% Pe: Process noise covariance matrix
% e:  residuals
% Pa: MVAR parameter covariance matrix
% y: (CH x LEN) data matrix
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

%% TODO:

% Replace 3rd dimension of Ah with Ah_. 
% Convert Ah_ to sparse and Ah
% Implement downsamping

%%

if verb
    h=waitbar(0,sprintf('fitting VAR[%d] model [mode=%s] ...', ...
        p, 'Dual EKF'));
end

%% 
CH = size(y,1);                            % Number of states (here, CH = N)
LEN = size(y,2);                           % Number of the multivariate observations

if nargin<3
    lambda = .02;                          % update coefficent for parameter covariance matrix (forgetting factor)
end

%% Initial parameters for Dual Extended Kalman Filter

if nargout>3
    err = zeros(size(y));
end

%%%%% (EKF 1)
xh = zeros(CH*p,LEN);                      % (EKF 1) Initial a-posteriori states (Mp x 1)
Px = .1*eye(CH*p);                         % (EKF 1) Initial a-posteriori state covariance matrix
R = eye(CH);
B = zeros(CH*p,CH);                        % (EKF 1) Relationship between states and the process noise ( x(k) = F[x(k-1)] + B*v(k) )
B(1:CH,:) = eye(CH);                       % (EKF 1) B = [I 0 ... 0]'

%%%%% EKF 2
Pa = eye(CH*CH*p);                         % (EKF 2) Initial a-posteriori parameters covariance matrix
Ah = zeros(CH*p,CH*p,LEN);                 % (EKF 2) Initial a-posteriori parameters estimates (matrix form of 'ah' plus identity matrices)
Ah(1:CH,1:CH*p,p) = .1*randn(CH,CH*p);     % (EKF 2) Initial a-posteriori parameters estimates (matrix form of 'ah' plus identity matrices)
for r = 2 : p
    Ah((r-1)*CH+1:r*CH,(r-2)*CH+1:(r-1)*CH, p) = eye(CH);
end

for r = 1 : p
    xh((r-1)*CH+1:r*CH,p+1) = y(:,p-r+1);
end

%%%%% Mutual variables between EKF 1 and DEKF 2
Q = 10*eye(CH);                            % (EKF 1,2) Initial process noise covariance matrix
C = B';                                    % (EKF 1,2) Measurement matrix (identity matrix, C = B')

Pe = zeros(CH,CH,LEN);                       % (EKF 1) Noise covariance matrix (TM)
Pe(:,:,1:p) = repmat(Q,[1,1,p]);

I_CHp       = eye(CH*p);
I_CH2p      = eye(CH*CH*p);
I_CHp_m1    = eye(CH*(p-1));

%% DEKF starts ....
for i = p+1 : LEN
      
    if verb && ~mod(i,floor(LEN/100))
        waitbar(i/LEN,h,...
             sprintf('fitting VAR[%d] model [mode=%s] (%d/%d) ...',p,'Dual EKF',i,LEN));
    end
    
    if any(isnan(y(:,i)))
        % don't update on missing data
        continue;
    end
      
    [J_x J_A] = MVAR_JacCSD(Ah(:,:,i-1),xh(:,i-1),p);  % xh(k) = F(A(k-1) * xh(k-1)) = Ah(k-1) * xh(k-1)
    Ah_ = Ah(:,:,i-1);                                 % Ah_(k) = Ah(k-1)

    %% EKF 1 ---> States estimation    
    %---------- Time Update (EKF1) ----------
    Rv = B * Q * B';                          % According to Haykin's book
    xh_ = Ah_ * xh(:,i-1);                    % xh_(k) = A_h(k-1) * xh(k-1)
    Px_ = J_x * Px * J_x' + Rv;               % Px_(k) = A_h(k-1) * Px(k-1) * A_h(k-1)' + B * Q * B'
    
    %---------- Measurement Update (EKF1) ----------
    Rn = R; %R(:,:,i-1);                      % According to Haykin's book
    Kx = Px_ * C' /(C * Px_ * C' + Rn);       % Kx(k)  = Px_(k) * C' * inv(C * Px_(k) * C' + R)
    Px = (I_CHp - Kx * C) * Px_;          % Px(k)  = (I - Kx(k) * C) * Px_(k)
    e = y(:,i) - C * xh_;                     % inov(k) = y(k) - C * Ah_(k) * xh(k-1)
    Pe(:,:,i) = (1-lambda)*Q + lambda*(e*e');                         % noise covariance matrix
    xh(:,i) = xh_ + Kx * e;                   % xh(k)  = xh_(k) + Kx(k) * (y(k) - C * xh_(k)) 
    
    if nargout>3
        err(:,i) = e;
    end
    
    %% EKF 2 ---> Parameters estimation
    %---------- Time Update (EKF2) ----------
    ah_ = reshape(Ah_(1:CH,:)',CH*CH*p,1);    % ah_ = vec(Ah(k-1))
    Rr = lambda*Pa;                           % Rr = lambda * Pa(k-1)
    Pa_ = Pa + Rr;                            % Pa_(k) = Pa(k-1) + Rr
    
    %---------- Measurement Update (EKF2) ----------
    %%%%%% Compute DfDa and H (H(k) = C*DfDa(k-1))
    H = C * (-J_A);                           % J_A = -DfDa;
    
    %%%%%%%%%%%%%%%
    Re = (Rn + Q);                            % According to Haykin's book
    Ka = Pa_ * H' /(H * Pa_ * H' + Re);       % Ka(k) = Pa_(k) * H(k) * inv(H(k) * Pa_(k) * H(k)' +R + Q)
    Pa = (I_CH2p - Ka * H) * Pa_;           % Pa(k) = (I - Ka(k) * H(k)) * Pa_(k)
    ah = ah_ + Ka * (y(:,i)- C * Ah_ * xh(:,i-1));    % ah(k) = ah_(k) + Ka(k) * (y(k) - yh_(k))
      
    %---------- Re-arrange vector ah(k) into the matrix Ah(k) ----------
    Ah(1:CH,1:CH*p,i) = reshape(ah,CH*p,CH)';
    Ah(CH+1:CH*p,1:(p-1)*CH,i) = I_CHp_m1;
%     for r = 2 : p
%         Ah((r-1)*CH+1:r*CH,(r-2)*CH+1:(r-1)*CH, i) = eye(CH);
%     end

end

A = Ah(1:CH,:,:); % Estimated time-varying parameters, A = [A1 A2 ... Ar]

if verb, close(h); end
