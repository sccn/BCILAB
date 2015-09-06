% V = est_calcInvCovMatFourierPDC(Rinv,E,foi,fs,N,p, verb)
%
% Obtain the frequency domain transform of the inverse covariance matrix of 
% an M-variate VAR[p] process
% This is needed for calculating the analytic statistics for the PDC [1]
%
% Inputs:
%
%     Rinv: inverse process covariance matrix obtained from est_calcInvCovMat
%     E:    noise covariance matrix
%     foi:  frequencies of interest (Hz)
%     fs:   sampling rate
%     N:    # chans
%     p:    model order
%     verb: verbosity level. 0 = no output, 1=text.
%
% Outputs:
%
%     V:    Frequency-domain transform of the inverse covariance matrix
%
% References:
%
% [1] Schelter et al, (2009). Testing for directed influences among neural 
% signals using partial directed coherence. J. Neuroscience Methods. 152:210-9.
%
% see also: est_calcInvCovMat(), est_calcInvCovMatFourier()
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function V = est_calcInvCovMatFourierPDC(Rinv,E,foi,fs,N,p,verb,DEBUG)

if nargin<8
    DEBUG = false;
end

%% extract the diagonal elements of H=Rinv
% structure of Hd is (e.g., for p=2): 
% [diag(H(1,1)), diag(H(2,1)), diag(H(1,2)), diag(H(2,2))]
% where diag(H(u,v)) is the column vector formed by the diagonal of
% submatrix H(u,v) of H
Hd = zeros(N,p^2);
cnt=1;
for v=1:p
    for u=1:p
        Hd(:,cnt)=diag(Rinv((u-1)*N+1:u*N,(v-1)*N+1:v*N));
        cnt=cnt+1;
    end
end

if DEBUG
    for j=1:N
        try chol(reshape(Hd(j,:),[p p])); catch; fprintf('Hj<0: j=%d - ',j); keyboard; end;
    end
end

%% create u,v index vectors
% us = [1 2 ... p 1 2 ... p ... ]'   (p^2 length)
% vs = [1 1 ... 1 2 2 ... 2 ... ]'   (p^2 length)
us = repmat(1:p,1,p)';
vs = zeros(p^2,1);
for ii=1:p
    vs((ii-1)*p+1:(ii-1)*p+p)=ones(p,1)*ii;
end



%% construct the V matrix for all freqs
fi=0;
COS=zeros(p^2,1); 
freqs=(2*pi*foi)/fs;
V = zeros(length(freqs),N,N); % NOTE: V will end up (N,N,freqs)
if verb, h=waitbar(0,'calculating V^-1...'); end
for f=freqs
    fi=fi+1;

    %% construct cosine matrix
    COS(:,1) = cos(us*f).*cos(vs*f)+sin(us*f).*sin(vs*f);
%     COS(:,2) = sin(us*f).*sin(vs*f);
%     COS(:,3) = cos(us*f).*sin(vs*f);
%     COS(:,4) = sin(us*f).*sin(vs*f);

    %% multiply Hd and COS matrices to get matrix where row j is 
    % [sum_{u,v=1 : p} Hjj(u,v)COS_11(u,v), ...
    %  sum_{u,v=1 : p} Hjj(u,v)COS_21(u,v), ...
    %  sum_{u,v=1 : p} Hjj(u,v)COS_12(u,v), ...
    %  sum_{u,v=1 : p} Hjj(u,v)COS_22(u,v)]
    %
    % where COSab(u,v) is the a,bth element of the sine-cos transform
    % matrix evaluated at u,v.
    %
    % NOTE: reshaping row j to 2x2 yeilds a matrix proportional to V_ij(f) 
    % for some specific i
    Vm = Hd*COS;
    
    % now we multiply in the variances of the i's (E(i,i)) to generate the
    % full V matrix
    Vm = kron(diag(E)',Vm);
    
    % next, reshape to desired structure (NxNx2x2)
    V(fi,:,:) = Vm';
    
    
    if verb, waitbar(fi/length(freqs),h); end
end % for freqs

% permute V: (chs,chs,freqs,2,2)
V = permute(V,[2 3 1]);
if verb, close(h); end

