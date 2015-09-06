function [xout,Q2out,e,Kalman, Kout] = mvaar(y,p,UC,mode,Kalman,verb,downsampleFactor,constraints)
% Multivariate (Vector) adaptive AR estimation base on a multidimensional
% Kalman filer algorithm. A standard VAR model (A0=I) is implemented. The
% state vector is defined as X=(A1|A2...|Ap)' and x=vec(X)
%
% [x,e,Kalman,Q2] = mvaar(y,p,UC,mode,Kalman)
%
% The standard MVAR model is defined as:
%
%		y(n)-A1(n)*y(n-1)-...-Ap(n)*y(n-p)=e(n)
%
%	The dimension of y(n) equals M
%
%	Input Parameters:
%
% 		y			Observed data or signal [N x M] where N = epoch length
% 		p			prescribed maximum model order (default 1)
%		UC			update coefficient	(default 0.001)
%		mode	 	update method of the process noise covariance matrix 0...7 ^
%					correspond to S0...S7 (default 0)
%       verb        verbosity
%       downsampleFactor:   Starting from sample t=max(2,k), store only every k Kalman
%                           coefficients (states, etc) where k=downsampleFactor.
%       constraints structure with fields .D and .d containing constraints
%                   of the form Dx = d (see [1] below)
%
%	Output Parameters
%
%		e			prediction error of dimension s
%		x			state matrix of dimension [T x M*M*p]
%                   where T = ceil((N-downsampleFactor+q)/downsampleFactor)
%                   where q = (downsampleFactor>1 ? 1 : 0)
%                   - note that we never store the coefficient matrix for t=1
%                   since this is always zero (we have no sample at t=0 from
%                   which to compute the coefficient matrix for t=1)
%		Q2			measurement noise covariance matrix of dimension M x M
%       Kout        estimated state noise covariance matrix
%       Kalman      Kalman structure (can be used as subsequent startup)
%

%       $Id: mvaar.m 5090 2008-06-05 08:12:04Z schloegl $
%       Copyright (C) 2001-2002 Christian Kasess
%       Copyright (C) 2003, 2008 Alois Schloegl
%
%       Copyright (C) 2010-2011 Tim Mullen
%
%       Modified by Tim Mullen
%       01/23/2011 -- Modified for downsampled storage
%       04/13/2011 -- Optimized performance
%       05/12/2011 -- Added projection onto constraint surface [1]
%       05/20/2011 -- Added additional noise covariance update modes
%
%       [1] Simon D (2010) Kalman Filtering with State Constraints:
%       A Survey of Linear and Nonlinear Algorithms.
%       Control Theory & Applications, IET 4:1303-€“1318
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

if nargin<6
    verb = 2;
end
if nargin<4,
    mode=0;
end;
if nargin<3,
    UC=0.001;
end;
if nargin<2,
    p=1;
end
if nargin<1,
    fprintf(2,'No arguments supplied\n');
    return
end;
if nargin<8 || isempty(constraints)
    doConstraints = 0;
else
    doConstraints = 1;
    Constr_D = constraints.D;
    Constr_d = constraints.d;
end

if ~any(mode==(0:8))
    fprintf(2,'Invalid mode (0...8)\n');
    return
end;


[LEN, M, NTR] = size(y);  % signal length, number of channels, number of trials
L = M*M*p;

if LEN<(p+1),
    fprintf(2,'Not enough observed data supplied for given model order\n');
    return
end

% ye = zeros(size(y));	%prediction of y

% size of downsampled storage
dslen = ceil((LEN-downsampleFactor+(downsampleFactor>1))/downsampleFactor);

xout=zeros(L,dslen,NTR);

if nargout>1
    Q2out=zeros(M,M,dslen,NTR);
end
if nargout>4
    Kout = zeros(M,M,dslen,NTR);
end;

if verb==2
    h=waitbar(0,sprintf('fitting VAR[%d] model [mode=%s] ...', ...
        p, 'Kalman'));
end

%Kalman Filter initialization (Kp (K predicted or a-priori) equals K(n+1,n) )
F   = eye(L);          % observation matrix
G   = zeros(L,M);      % Kalman Gain
x   = zeros(L,1);      % state vector
Kp  = eye(L);
Q1  = eye(L);         % state noise covariance matrix
Q2  = eye(M);         % measurement noise covariance matrix
ye  = zeros(size(y)); % prediction of y

if nargin>=4 && ~isempty(Kalman)
    
    if isfield(Kalman,'ye'), ye = Kalman.ye; end
    if isfield(Kalman,'F'),  F  = Kalman.F;  end
    if isfield(Kalman,'Q1'), Q1 = Kalman.Q1; end
    if isfield(Kalman,'Kp'), Kp = Kalman.Kp; end
    if isfield(Kalman,'Q2'), Q2 = Kalman.Q2; end
    if isfield(Kalman,'x'),  x  = Kalman.x;  end
    if isfield(Kalman,'H'),  H  = Kalman.H;  end
    if isfield(Kalman,'G'),  G  = Kalman.G;  end
    
end

upd = eye(L)/L*UC;		%diagonal matrix containing UC

if(mode==3)
    Block=kron(eye(M),ones(M*p));
elseif(mode==4)
    index=[];
    Block1=[];
    Block0=[];
    for i=1:M,
        index=[index ((i-1)*M*p+i:M:i*M*p)];
        mone=eye(M);
        mone(i,i)=0;
        mzero=eye(M)-mone;
        Block1=blkdiag(Block1,kron(eye(p),mone));
        Block0=blkdiag(Block0,kron(eye(p),mzero));
    end;
elseif mode==5
    Q1 = upd;  % a4 of thesis
elseif mode==6
    Q1 = eye(L)*UC^2;
elseif mode==7
    Q1 = eye(L)*UC;
elseif mode==8
    Q1 = zeros(L);  % RLS algorithm (no process noise)
end



for tr=1:NTR
    
    % TODO: add option to re-initialize the state variables here
    
    curval = 1;
    
    for n = 2:LEN
        
        if verb==2 && ~mod(n,100)
            waitbar(n/LEN,h,...
                {sprintf('Trial (%d/%d)',tr,NTR), ...
                 sprintf('fitting VAR[%d] model [mode=%s] (%d/%d) ...',p,'Kalman',n,LEN)});
        end
        
        
        if(n<=p)
            Yr=[y(n-1:-1:1,:,tr)' zeros(M,p-n+1)];	%vector of past observations
            Yr=Yr(:)';
        else
            Yr=y(n-1:-1:n-p,:,tr)';				%vector of past observations
            Yr=Yr(:)';
        end
        
        %Update of measurement matrix
        H = blkdiageye(Yr,M);        
        
        %calculate a-priori prediction error (1-step prediction error)
        ye(n,:,tr)=(H*x)';
        err=y(n,:,tr)-ye(n,:,tr);
        
        if ~any(isnan(err(:))),
            % update of Q2 (measurement noise covariance matrix, V)) 
            % using the prediction error of the previous step
            Q2=(1-UC)*Q2+UC*(err'*err);
            
            
            KpH=Kp*H';
            HKp=H*Kp;
            
            % Kalman gain
            G=KpH/(H*KpH+Q2);
            
            % calculation of the a-posteriori state error covariance matrix
            % K=Kp-G*KpH'; Althouh PK is supposed to be symmetric, this operation makes the filter unstable
            K=Kp-G*HKp;
            
            % mode==0 no update of Q1 (process noise covariance matrix, W)
            % update of Q1 using the predicted state error cov matrix
            if (mode==1)
                Q1=diag(diag(K)).*UC;
            elseif(mode==2)  % similar to mode a2 from thesis
                Q1=upd*trace(K);
            elseif(mode==3)
                Q1=diag(sum((Block*diag(diag(K))),2)')/(p*M)*UC;
            elseif(mode==4)
                avg=trace(K(index,index))/(p*M)*UC;
                Q1=Block1*UC+Block0*avg;
            end
            
            %a-priori state error covariance matrix for the next time step
            Kp=K+Q1;
            
            %current estimation of state x
            x=x+G*(err)';
            
            % perform optional constraint projection
            if doConstraints
                KD = K*Constr_D';
                
                % project the solution onto the constraint surface
                x = x - (KD/(Constr_D*KD))*(Constr_D*x - Constr_d);
            end
            
        end; % isnan(err)
        
        if ~mod(n,downsampleFactor)
            
            % store the current state
            xout(:,curval,tr) = x;
            
            if nargout>1
                Q2out(:,:,curval,tr)=Q2;
            end
            
            if nargout>4
                Kout(:,:,curval,tr) = K;
            end
            
            curval = curval + 1;
        end
    end;
    
end

if nargout > 2
    e = y - ye;
end

if nargout > 3
    Kalman.ye = ye;
    Kalman.F  = F;
    Kalman.Q1 = Q1;
    Kalman.Kp = Kp;
    Kalman.Q2 = Q2;
    Kalman.x = x;
    Kalman.H = H;
    Kalman.G = G;
end

xout = permute(xout,[2 1 3]);

if verb==2, close(h); end